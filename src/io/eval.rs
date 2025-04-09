use super::{params::ObfuscationParams, utils::*, Obfuscation};
use crate::{
    bgg::{sampler::BGGPublicKeySampler, BggEncoding, BitToInt},
    poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams},
};
use itertools::Itertools;
use std::sync::Arc;

impl<M> Obfuscation<M>
where
    M: PolyMatrix,
{
    pub fn eval<SH, ST>(
        &self,
        obf_params: ObfuscationParams<M>,
        mut sampler_hash: SH,
        inputs: &[bool],
    ) -> Vec<bool>
    where
        SH: PolyHashSampler<[u8; 32], M = M>,
        ST: PolyTrapdoorSampler<M = M>,
    {
        sampler_hash.set_key(self.hash_key);
        let params = Arc::new(obf_params.params.clone());
        let d = obf_params.d;
        let d1 = d + 1;
        let sampler = Arc::new(sampler_hash);
        debug_assert_eq!(inputs.len(), obf_params.input_size);
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone(), d);
        let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
        let packed_input_size = public_data.packed_input_size;
        let packed_output_size = public_data.packed_output_size;
        let (mut ps, mut encodings) = (vec![], vec![]);
        ps.push(self.p_init.clone());
        encodings.push(self.encodings_init.clone());

        // Sample public keys
        #[cfg(feature = "test")]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
        #[cfg(not(feature = "test"))]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();
        let pubkeys = (0..obf_params.input_size + 1)
            .map(|idx| {
                sample_public_key_by_idx(
                    &bgg_pubkey_sampler,
                    &obf_params.params,
                    idx,
                    &reveal_plaintexts,
                )
            })
            .collect_vec();

        #[cfg(feature = "test")]
        if obf_params.encoding_sigma == 0.0 &&
            obf_params.hardcoded_key_sigma == 0.0 &&
            obf_params.p_sigma == 0.0
        {
            let expected_p_init = {
                let s_connect = self.s_init.concat_columns(&[&self.s_init]);
                s_connect * &self.bs[0][2]
            };
            debug_assert_eq!(self.p_init, expected_p_init);

            let zero = <M::P as Poly>::const_zero(&params);
            let one = <M::P as Poly>::const_one(&params);
            let inserted_poly_gadget = {
                let mut polys = vec![];
                polys.push(one.clone());
                for _ in 0..(obf_params.input_size.div_ceil(params.ring_dimension() as usize)) {
                    polys.push(zero.clone());
                }
                polys.push(self.t_bar.clone());
                let gadget_d1 = M::gadget_matrix(&params, d1);
                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_d1)
            };
            let expected_encoding_init = self.s_init.clone() *
                &(pubkeys[0][0].concat_matrix(&pubkeys[0][1..]) - inserted_poly_gadget);
            debug_assert_eq!(
                encodings[0][0].concat_vector(&encodings[0][1..]),
                expected_encoding_init
            );
        }
        let log_base_q = params.as_ref().modulus_digits();
        let dim = params.as_ref().ring_dimension() as usize;
        for (idx, input) in inputs.iter().enumerate() {
            let m = if *input { &self.m_preimages[idx][1] } else { &self.m_preimages[idx][0] };
            let q = ps[idx].clone() * m;
            let n = if *input { &self.n_preimages[idx][1] } else { &self.n_preimages[idx][0] };
            let p = q.clone() * n;
            let k = if *input { &self.k_preimages[idx][1] } else { &self.k_preimages[idx][0] };
            let v = q.clone() * k;
            let new_encode_vec = {
                let t = if *input { &public_data.rgs[1] } else { &public_data.rgs[0] };
                let encode_vec = encodings[idx][0].concat_vector(&encodings[idx][1..]);
                let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
                encode_vec.mul_tensor_identity_decompose(t, packed_input_size + 1) + v
            };
            let mut new_encodings = vec![];
            let inserted_poly_index = 1 + idx / dim;
            for (j, encode) in encodings[idx].iter().enumerate() {
                let m = d1 * log_base_q;
                let new_vec = new_encode_vec.slice_columns(j * m, (j + 1) * m);
                let plaintext = if j == inserted_poly_index {
                    let inserted_coeff_index = idx % dim;
                    let mut coeffs = encode.plaintext.as_ref().unwrap().coeffs().clone();
                    coeffs[inserted_coeff_index] = if *input {
                        <M::P as Poly>::Elem::one(&params.modulus())
                    } else {
                        <M::P as Poly>::Elem::zero(&params.modulus())
                    };
                    Some(M::P::from_coeffs(params.as_ref(), &coeffs))
                } else {
                    encode.plaintext.clone()
                };
                let new_pubkey = pubkeys[idx + 1][j].clone();
                let new_encode: BggEncoding<M> =
                    BggEncoding::new(new_vec, new_pubkey.clone(), plaintext);
                new_encodings.push(new_encode);
            }
            ps.push(p.clone());
            encodings.push(new_encodings);
            #[cfg(feature = "test")]
            if obf_params.encoding_sigma == 0.0 &&
                obf_params.hardcoded_key_sigma == 0.0 &&
                obf_params.p_sigma == 0.0
            {
                let mut cur_s = self.s_init.clone();
                for bit in inputs[0..idx].iter() {
                    let r = if *bit { public_data.r_1.clone() } else { public_data.r_0.clone() };
                    cur_s = cur_s * r;
                }
                let new_s = if *input {
                    cur_s.clone() * &public_data.r_1
                } else {
                    cur_s.clone() * &public_data.r_0
                };
                let b_next_bit =
                    if *input { self.bs[idx + 1][1].clone() } else { self.bs[idx + 1][0].clone() };
                let expected_q = cur_s.concat_columns(&[&new_s]) * &b_next_bit;
                debug_assert_eq!(q, expected_q);
                let expected_p = new_s.concat_columns(&[&new_s]) * &self.bs[idx + 1][2];
                debug_assert_eq!(p, expected_p);
                let expcted_new_encode = {
                    let dim = params.ring_dimension() as usize;
                    let one = <M::P as Poly>::const_one(&params);
                    let gadget_d1 = M::gadget_matrix(&params, d1);
                    let inserted_poly_gadget = {
                        let mut polys = vec![];
                        polys.push(one.clone());
                        let mut coeffs = vec![];
                        for bit in inputs[0..=idx].iter() {
                            if *bit {
                                coeffs.push(<M::P as Poly>::Elem::one(&params.modulus()));
                            } else {
                                coeffs.push(<M::P as Poly>::Elem::zero(&params.modulus()));
                            }
                        }
                        for _ in 0..(obf_params.input_size - idx - 1) {
                            coeffs.push(<M::P as Poly>::Elem::zero(&params.modulus()));
                        }
                        let input_polys = coeffs
                            .chunks(dim)
                            .map(|coeffs| M::P::from_coeffs(&params, coeffs))
                            .collect_vec();
                        polys.extend(input_polys);
                        polys.push(self.t_bar.clone());
                        M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_d1)
                    };
                    let pubkey = pubkeys[idx + 1][0].concat_matrix(&pubkeys[idx + 1][1..]);
                    new_s * (pubkey - inserted_poly_gadget)
                };
                debug_assert_eq!(new_encode_vec, expcted_new_encode);
            }
        }
        let a_decomposed = public_data.a_rlwe_bar.entry(0, 0).decompose_bits(params.as_ref());
        let b_decomposed = &self.ct_b.entry(0, 0).decompose_bits(params.as_ref());
        let final_circuit = build_final_bits_circuit::<M::P, BggEncoding<M>>(
            &a_decomposed,
            b_decomposed,
            obf_params.public_circuit.clone(),
        );
        let last_input_encodings = encodings.last().unwrap();
        let output_encodings = final_circuit.eval::<BggEncoding<M>>(
            params.as_ref(),
            &last_input_encodings[0],
            &last_input_encodings[1..],
        );
        let log_q = params.as_ref().modulus_bits();
        let output_encoding_ints = output_encodings
            .chunks(log_q)
            .map(|bits| BggEncoding::bits_to_int(bits, &params))
            .collect_vec();
        let output_encodings_vec =
            output_encoding_ints[0].concat_vector(&output_encoding_ints[1..]);
        let final_preimage = &self.final_preimage;
        let final_v = ps.last().unwrap().clone() * final_preimage;
        let z = output_encodings_vec.clone() - final_v.clone();
        debug_assert_eq!(z.size(), (1, packed_output_size));
        #[cfg(feature = "test")]
        if obf_params.encoding_sigma == 0.0 &&
            obf_params.hardcoded_key_sigma == 0.0 &&
            obf_params.p_sigma == 0.0
        {
            let mut last_s = self.s_init.clone();
            for bit in inputs.iter() {
                let r = if *bit { public_data.r_1.clone() } else { public_data.r_0.clone() };
                last_s = last_s * r;
            }

            let output_plaintext =
                output_encoding_ints[0].plaintext.as_ref().unwrap().extract_highest_bits();
            let hardcoded_key_bits = self
                .hardcoded_key
                .coeffs()
                .iter()
                .map(|elem| elem != &<M::P as Poly>::Elem::zero(&params.modulus()))
                .collect::<Vec<_>>();
            debug_assert_eq!(output_plaintext, hardcoded_key_bits);
            {
                let expcted = last_s *
                    (output_encoding_ints[0].pubkey.matrix.clone() -
                        M::unit_column_vector(params.as_ref(), d1, d1 - 1) *
                            output_encoding_ints[0].plaintext.clone().unwrap());
                debug_assert_eq!(output_encoding_ints[0].vector, expcted);
            }
            debug_assert_eq!(z.size(), (1, packed_output_size));
            debug_assert_eq!(z.entry(0, 0), output_encoding_ints[0].plaintext.clone().unwrap());
        }
        z.get_row(0).into_iter().flat_map(|p| p.extract_bits_with_threshold(&params)).collect_vec()
    }
}
