use super::utils::*;
use super::{Obfuscation, ObfuscationParams};
use crate::bgg::sampler::BGGPublicKeySampler;
use crate::bgg::BggEncoding;
use crate::poly::{matrix::*, sampler::*, Poly, PolyParams};
use itertools::Itertools;
use std::sync::Arc;

pub fn eval_obf<M, SH>(
    obf_params: ObfuscationParams<M>,
    mut sampler_hash: SH,
    obfuscation: Obfuscation<M>,
    input: &[bool],
) -> Vec<bool>
where
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M>,
{
    sampler_hash.set_key(obfuscation.hash_key);
    let params = Arc::new(obf_params.params.clone());
    let sampler = Arc::new(sampler_hash);
    // let dim = params.as_ref().ring_dimension() as usize;
    debug_assert_eq!(input.len(), obf_params.input_size);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler.clone());
    let public_data = PublicSampledData::sample(&obf_params, &bgg_pubkey_sampler);
    // let packed_input_size = public_data.packed_input_size;
    let packed_output_size = public_data.packed_output_size;
    let (mut ps, mut encodings) = (vec![], vec![]);
    ps.push(obfuscation.p_init.clone());
    encodings.push(obfuscation.encodings_init);
    // let encode_inputs =
    //     obfuscation.encode_input.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
    // cs_input.push(encode_inputs[0].concat_columns(&encode_inputs[1..]));
    // let encode_fhe_key =
    //     obfuscation.encode_fhe_key.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
    // cs_fhe_key.push(encode_fhe_key[0].concat_columns(&encode_fhe_key[1..]));
    let log_q = params.as_ref().modulus_bits();
    for (idx, input) in input.iter().enumerate() {
        let m =
            if *input { &obfuscation.m_preimages[idx].1 } else { &obfuscation.m_preimages[idx].0 };
        let q = ps[idx].clone() * m;
        let n =
            if *input { &obfuscation.n_preimages[idx].1 } else { &obfuscation.n_preimages[idx].0 };
        let p = q.clone() * n;
        let k =
            if *input { &obfuscation.k_preimages[idx].1 } else { &obfuscation.k_preimages[idx].0 };
        let v = q * k;
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
            let new_encode = BggEncoding::new(new_vec, new_pubkey.clone(), plaintext);
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
        ps.push(p);
        encodings.push(new_encodings);
        // cs_input.push(c_input);
        // cs_fhe_key.push(c_fhe_key);
    }
    let final_circuit =
        build_final_step_circuit::<_, BggEncoding<M>>(&params, obf_params.public_circuit.clone());
    let last_input_encodings = encodings.last().unwrap();
    let output_encodings = final_circuit.eval_poly_circuit::<BggEncoding<M>>(
        &params,
        last_input_encodings[0].clone(),
        &last_input_encodings[1..],
    );
    let final_v = ps.last().unwrap().clone() * &obfuscation.final_preimage;
    let z = output_encodings[0].concat_vector(&output_encodings[1..]) - final_v;
    debug_assert_eq!(z.size(), (1, packed_output_size));
    z.get_row(0).into_iter().flat_map(|p| p.extract_highest_bits()).collect_vec()
}
