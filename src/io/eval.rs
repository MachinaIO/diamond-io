use super::{utils::*, Obfuscation, ObfuscationParams};
use crate::{
    bgg::{sampler::BGGPublicKeySampler, BggEncoding},
    poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams},
};
use itertools::Itertools;
use std::sync::Arc;

impl<M> Obfuscation<M>
where
    M: PolyMatrix,
{
    pub fn eval<SH>(
        &self,
        obf_params: ObfuscationParams<M>,
        mut sampler_hash: SH,
        inputs: &[bool],
    ) -> Vec<bool>
    where
        SH: PolyHashSampler<[u8; 32], M = M>,
    {
        sampler_hash.set_key(self.hash_key);
        let params = Arc::new(obf_params.params.clone());
        let sampler = Arc::new(sampler_hash);
        debug_assert_eq!(inputs.len(), obf_params.input_size);
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone());
        let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
        let packed_output_size = public_data.packed_output_size;
        let (mut ps, mut encodings) = (vec![], vec![]);
        ps.push(self.p_init.clone());
        encodings.push(self.encodings_init.clone());
        #[cfg(test)]
        {
            let expected_p_init = {
                let s_connect = self.s_init.clone().concat_columns(&[self.s_init.clone()]);
                s_connect * &self.bs[0].2
            };
            debug_assert_eq!(self.p_init, expected_p_init);

            let zero = <M::P as Poly>::const_zero(&params);
            let one = <M::P as Poly>::const_one(&params);
            let enc_hardcoded_key_decomposed = self.enc_hardcoded_key.decompose().get_column(0);
            let inserted_poly_gadget = {
                let mut polys = vec![];
                polys.push(one.clone());
                for poly in enc_hardcoded_key_decomposed.iter() {
                    polys.push(poly.clone());
                }
                for _ in 0..(obf_params.input_size.div_ceil(params.ring_dimension() as usize)) {
                    polys.push(zero.clone());
                }
                polys.push(self.t_bar.entry(0, 0).clone());
                let gadget_2 = M::gadget_matrix(&params, 2);
                M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_2)
            };
            let expected_encoding_init = self.s_init.clone() *
                (public_data.pubkeys[0][0].concat_matrix(&public_data.pubkeys[0][1..]) -
                    inserted_poly_gadget);
            debug_assert_eq!(
                encodings[0][0].concat_vector(&encodings[0][1..]),
                expected_encoding_init
            );
        }
        // let encode_inputs =
        //     self.encode_input.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
        // cs_input.push(encode_inputs[0].concat_columns(&encode_inputs[1..]));
        // let encode_fhe_key =
        //     self.encode_fhe_key.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
        // cs_fhe_key.push(encode_fhe_key[0].concat_columns(&encode_fhe_key[1..]));
        let log_q = params.as_ref().modulus_bits();
        let dim = params.as_ref().ring_dimension() as usize;
        for (idx, input) in inputs.iter().enumerate() {
            let m = if *input { &self.m_preimages[idx].1 } else { &self.m_preimages[idx].0 };
            let q = ps[idx].clone() * m;
            let n = if *input { &self.n_preimages[idx].1 } else { &self.n_preimages[idx].0 };
            let p = q.clone() * n;
            let k = if *input { &self.k_preimages[idx].1 } else { &self.k_preimages[idx].0 };
            let v = q.clone() * k;
            // let v_input = v.slice_columns(0, 2 * log_q * (packed_input_size + 1));
            // let v_fhe_key = v.slice_columns(
            //     2 * log_q * (packed_input_size + 1),
            //     2 * log_q * (packed_input_size + 3),
            // );
            let new_encode_vec = {
                let t = if *input {
                    &public_data.rgs_decomposed[1]
                } else {
                    &public_data.rgs_decomposed[0]
                };
                let encode_vec = encodings[idx][0].concat_vector(&encodings[idx][1..]);
                let packed_input_size = log_q + obf_params.input_size.div_ceil(dim) + 1;
                encode_vec.mul_tensor_identity(t, packed_input_size + 1) + v
            };
            let mut new_encodings = vec![];
            // let zero_poly = <M::P as Poly>::const_zero(&params);
            // let one_poly = <M::P as Poly>::const_one(&params);
            let inserted_poly_index = 1 + log_q + idx / dim;
            for (j, encode) in encodings[idx].iter().enumerate() {
                // let encode = encodings[idx][j].clone();
                let m = 2 * log_q;
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
                let new_pubkey = public_data.pubkeys[idx + 1][j].clone();
                let new_encode: BggEncoding<M> =
                    BggEncoding::new(new_vec, new_pubkey.clone(), plaintext);
                new_encodings.push(new_encode);
            }
            // let c_input = {
            //     let t = if *input { &public_data.t_1.0 } else { &public_data.t_0.0 };
            //     cs_input[idx].clone() * t + v_input
            // };
            // let c_fhe_key = {
            //     let t = if *input { &public_data.t_1.1 } else { &public_data.t_0.1 };
            //     cs_fhe_key[idx].clone() * t + v_fhe_key
            // };
            ps.push(p.clone());
            encodings.push(new_encodings);
            #[cfg(test)]
            {
                let mut cur_s = self.s_init.clone();
                for bit in inputs[0..idx].iter() {
                    let r = if *bit { public_data.r_1.clone() } else { public_data.r_0.clone() };
                    cur_s = cur_s * r;
                }
                let new_s = if *input {
                    cur_s.clone() * public_data.r_1.clone()
                } else {
                    cur_s.clone() * public_data.r_0.clone()
                };
                let b_next_bit =
                    if *input { self.bs[idx + 1].1.clone() } else { self.bs[idx + 1].0.clone() };
                let expected_q = cur_s.concat_columns(&[new_s.clone()]) * &b_next_bit;
                debug_assert_eq!(q, expected_q);
                let expected_p = new_s.concat_columns(&[new_s.clone()]) * &self.bs[idx + 1].2;
                debug_assert_eq!(p, expected_p);
                let expcted_new_encode = {
                    let dim = params.ring_dimension() as usize;
                    let one = <M::P as Poly>::const_one(&params);
                    let gadget_2 = M::gadget_matrix(&params, 2);
                    let enc_hardcoded_key_decomposed =
                        self.enc_hardcoded_key.decompose().get_column(0);
                    let inserted_poly_gadget = {
                        let mut polys = vec![];
                        polys.push(one.clone());
                        for poly in enc_hardcoded_key_decomposed.iter() {
                            polys.push(poly.clone());
                        }
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
                        polys.push(self.t_bar.entry(0, 0).clone());
                        M::from_poly_vec_row(params.as_ref(), polys).tensor(&gadget_2)
                    };
                    let pubkey = public_data.pubkeys[idx + 1][0]
                        .concat_matrix(&public_data.pubkeys[idx + 1][1..]);
                    new_s * (pubkey - inserted_poly_gadget)
                };
                debug_assert_eq!(new_encode_vec, expcted_new_encode);
            }
            // cs_input.push(c_input);
            // cs_fhe_key.push(c_fhe_key);
        }
        let a_decomposed_polys = public_data.a_rlwe_bar.decompose().get_column(0);
        let final_circuit = build_final_step_circuit::<_, BggEncoding<M>>(
            &params,
            &a_decomposed_polys,
            obf_params.public_circuit.clone(),
        );
        let last_input_encodings = encodings.last().unwrap();
        let output_encodings = final_circuit.eval::<BggEncoding<M>>(
            &params,
            last_input_encodings[0].clone(),
            &last_input_encodings[1..],
        );
        let identity_2 = M::identity(&params, 2, None);
        let unit_vector = identity_2.slice_columns(1, 2);
        let output_encodings_vec =
            output_encodings[0].concat_vector(&output_encodings[1..]) * unit_vector.decompose();
        let final_v = ps.last().unwrap().clone() * &self.final_preimage;
        let z = output_encodings_vec.clone() - final_v.clone();
        debug_assert_eq!(z.size(), (1, packed_output_size));
        #[cfg(test)]
        {
            let mut last_s = self.s_init.clone();
            for bit in inputs.iter() {
                let r = if *bit { public_data.r_1.clone() } else { public_data.r_0.clone() };
                last_s = last_s * r;
            }

            let output_plaintext =
                output_encodings[0].plaintext.as_ref().unwrap().extract_highest_bits();
            let hardcoded_key_bits = self
                .hardcoded_key
                .entry(0, 0)
                .coeffs()
                .iter()
                .map(|elem| elem != &<M::P as Poly>::Elem::zero(&params.modulus()))
                .collect::<Vec<_>>();
            debug_assert_eq!(output_plaintext, hardcoded_key_bits);
            {
                let expcted = last_s.clone() *
                    (output_encodings[0].pubkey.matrix.clone() -
                        M::gadget_matrix(&params, 2) *
                            output_encodings[0].plaintext.clone().unwrap());
                debug_assert_eq!(output_encodings[0].vector, expcted);
            }

            // let a_f = obfuscation.final_preimage_target.slice_rows(0, 2).clone();
            // debug_assert_eq!(output_encodings[0].pubkey.matrix.clone() * unit_vector.decompose(),
            // a_f); let expected_final_v = last_s.clone() * &a_f;
            // debug_assert_eq!(final_v, expected_final_v);

            // let expected_output_encodings_vec = last_s.clone() * &a_f
            //     + M::from_poly_vec_row( &params,
            //       vec![output_encodings[0].plaintext.as_ref().unwrap().clone()],
            //     );
            // debug_assert_eq!(output_encodings_vec, expected_output_encodings_vec);

            // let scale = M::P::from_const(&params, &<M::P as
            // Poly>::Elem::half_q(&params.modulus())); debug_assert_eq!(z,
            // obfuscation.hardcoded_key * scale);
        }
        z.get_row(0).into_iter().flat_map(|p| p.extract_highest_bits()).collect_vec()
    }
}
