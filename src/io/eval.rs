use super::utils::*;
use super::{Obfuscation, ObfuscationParams};
use crate::bgg::sampler::BGGPublicKeySampler;
use crate::bgg::BggEncoding;
use crate::poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams};
use itertools::Itertools;
use std::sync::Arc;

pub fn eval_obf<M, SH>(
    obf_params: ObfuscationParams<M>,
    mut sampler_hash: SH,
    obfuscation: Obfuscation<M>,
    inputs: &[bool],
) -> Vec<bool>
where
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M>,
{
    sampler_hash.set_key(obfuscation.hash_key);
    let params = Arc::new(obf_params.params.clone());
    let sampler = Arc::new(sampler_hash);
    // let dim = params.as_ref().ring_dimension() as usize;
    debug_assert_eq!(inputs.len(), obf_params.input_size);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone());
    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    // let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let (mut ps, mut encodings) = (vec![], vec![]);
    ps.push(obfuscation.p_init.clone());
    encodings.push(obfuscation.encodings_init);
    #[cfg(test)]
    {
        let expected_p_init = {
            let s_connect =
                obfuscation.s_init.clone().concat_columns(&[obfuscation.s_init.clone()]);
            s_connect * &obfuscation.bs[0].2
        };
        debug_assert_eq!(obfuscation.p_init, expected_p_init);
    }
    // let encode_inputs =
    //     obfuscation.encode_input.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
    // cs_input.push(encode_inputs[0].concat_columns(&encode_inputs[1..]));
    // let encode_fhe_key =
    //     obfuscation.encode_fhe_key.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
    // cs_fhe_key.push(encode_fhe_key[0].concat_columns(&encode_fhe_key[1..]));
    let log_q = params.as_ref().modulus_bits();
    for (idx, input) in inputs.iter().enumerate() {
        let m =
            if *input { &obfuscation.m_preimages[idx].1 } else { &obfuscation.m_preimages[idx].0 };
        let q = ps[idx].clone() * m;
        let n =
            if *input { &obfuscation.n_preimages[idx].1 } else { &obfuscation.n_preimages[idx].0 };
        let p = q.clone() * n;
        let k =
            if *input { &obfuscation.k_preimages[idx].1 } else { &obfuscation.k_preimages[idx].0 };
        let v = q.clone() * k;
        // let v_input = v.slice_columns(0, 2 * log_q * (packed_input_size + 1));
        // let v_fhe_key = v.slice_columns(
        //     2 * log_q * (packed_input_size + 1),
        //     2 * log_q * (packed_input_size + 3),
        // );
        let new_encode_vec = {
            let t = if *input { &public_data.ts[1] } else { &public_data.ts[0] };
            let encode_vec = encodings[idx][0].concat_vector(&encodings[idx][1..]);
            encode_vec * t + v
        };
        let mut new_encodings = vec![];
        let zero_poly = <M::P as Poly>::const_zero(&params);
        let one_poly = <M::P as Poly>::const_one(&params);
        for (j, new_pubkey) in public_data.pubkeys[idx].iter().enumerate() {
            let encode = encodings[idx][j].clone();
            let m = 2 * log_q;
            let new_vec = new_encode_vec.slice_columns(j * m, (j + 1) * m);
            let plaintext = if j == 1 + log_q + idx {
                if *input {
                    Some(one_poly.clone())
                } else {
                    Some(zero_poly.clone())
                }
            } else {
                encode.plaintext.clone()
            };
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
            let mut cur_s = obfuscation.s_init.clone();
            for j in 0..idx {
                let r = if inputs[j] { public_data.r_1.clone() } else { public_data.r_0.clone() };
                cur_s = cur_s * r;
            }
            let new_s = if *input {
                cur_s.clone() * public_data.r_1.clone()
            } else {
                cur_s.clone() * public_data.r_0.clone()
            };
            let b_next_bit = if *input {
                obfuscation.bs[idx + 1].1.clone()
            } else {
                obfuscation.bs[idx + 1].0.clone()
            };
            let expected_q = cur_s.concat_columns(&[new_s.clone()]) * &b_next_bit;
            debug_assert_eq!(q, expected_q);
            let expected_p = new_s.concat_columns(&[new_s.clone()]) * &obfuscation.bs[idx + 1].2;
            debug_assert_eq!(p, expected_p);
            let expcted_new_encode = {
                let dim = params.ring_dimension() as usize;
                // let inserted_poly_index = 1 + log_q + idx / dim;
                // let inserted_coeff_index = idx % dim;
                // let zero_coeff = <M::P as Poly>::Elem::zero(&params.modulus());
                // let mut coeffs = vec![zero_coeff; dim];
                // coeffs[inserted_coeff_index] = if !input {
                //     <M::P as Poly>::Elem::zero(&params.modulus())
                // } else {
                //     <M::P as Poly>::Elem::one(&params.modulus())
                // };
                // let zero = <M::P as Poly>::const_zero(&params);
                let one = <M::P as Poly>::const_one(&params);
                // let inserted_poly = M::P::from_coeffs(params.as_ref(), &coeffs);
                let gadget_2 = M::gadget_matrix(&params, 2);
                let a_rlwe_decomposed = public_data.a_rlwe_bar.decompose().get_column(0);
                let inserted_poly_gadget = {
                    let mut polys = vec![];
                    polys.push(one.clone());
                    for j in 0..log_q {
                        polys.push(a_rlwe_decomposed[j].clone());
                    }
                    let mut coeffs = vec![];
                    for j in 0..=idx {
                        if inputs[j] {
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
                    println!("input_polys {:?}", input_polys);
                    polys.extend(input_polys);
                    polys.push(obfuscation.t_bar.entry(0, 0).clone());
                    // for _ in 0..(inserted_poly_index) {
                    //     polys.push(zero.clone());
                    // }
                    // polys.push(inserted_poly);
                    // for _ in (inserted_poly_index + 1)..(public_data.packed_input_size + 1) {
                    //     polys.push(zero.clone());
                    // }
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
    {
        let public_circuit_output = obf_params.public_circuit.clone().eval(
            &params,
            last_input_encodings[0].clone(),
            &last_input_encodings[1..last_input_encodings.len() - 1],
        );
        println!("public_circuit_output {:?}", public_circuit_output);
    }
    let output_encodings = final_circuit.eval::<BggEncoding<M>>(
        &params,
        last_input_encodings[0].clone(),
        &last_input_encodings[1..],
    );
    // println!("output_encodings {:?}", output_encodings);
    let identity_2 = M::identity(&params, 2, None);
    let unit_vector = identity_2.slice_columns(1, 2);
    let output_encodings_vec =
        output_encodings[0].concat_vector(&output_encodings[1..]) * unit_vector.decompose();
    let final_v = ps.last().unwrap().clone() * &obfuscation.final_preimage;
    let z = output_encodings_vec - final_v;
    debug_assert_eq!(z.size(), (1, packed_output_size));
    z.get_row(0).into_iter().flat_map(|p| p.extract_highest_bits()).collect_vec()
}
