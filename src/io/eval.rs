use super::utils::*;
use super::Obfuscation;
use crate::bgg::sampler::BGGPublicKeySampler;
use crate::poly::{matrix::*, sampler::*, Poly, PolyParams};
use crate::utils::*;
use itertools::Itertools;
use std::sync::Arc;

pub fn eval_obf<M, S>(
    params: <M::P as Poly>::Params,
    mut sampler: S,
    obfuscation: Obfuscation<M>,
    input: &[bool],
) -> Vec<bool>
where
    M: PolyMatrix,
    S: PolyHashSampler<[u8; 32], M = M>,
{
    sampler.set_key(obfuscation.hash_key);
    let params = Arc::new(params);
    let sampler = Arc::new(sampler);
    let dim = params.as_ref().ring_dimension() as usize;
    let input_size = input.len();
    let packed_input_size = input_size.div_ceil(dim);
    let bgg_pubkey_sampler = BGGPublicKeySampler::new(params.clone(), sampler.clone());
    let public_data = PublicSampledData::sample(
        params.as_ref(),
        sampler,
        &bgg_pubkey_sampler,
        input_size,
        packed_input_size,
        1,
    );
    let (mut ps, mut cs_input, mut cs_fhe_key) = (vec![], vec![], vec![]);
    ps.push(obfuscation.p_init.clone());
    let encode_inputs =
        obfuscation.encode_input.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
    cs_input.push(encode_inputs[0].concat_columns(&encode_inputs[1..]));
    let encode_fhe_key =
        obfuscation.encode_fhe_key.iter().map(|pubkey| pubkey.vector.clone()).collect_vec();
    cs_fhe_key.push(encode_fhe_key[0].concat_columns(&encode_fhe_key[1..]));
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
        let v_input = v.slice_columns(0, 2 * log_q * (packed_input_size + 1));
        let v_fhe_key = v.slice_columns(
            2 * log_q * (packed_input_size + 1),
            2 * log_q * (packed_input_size + 3),
        );
        let c_input = {
            let t = if *input { &public_data.t_1.0 } else { &public_data.t_0.0 };
            cs_input[idx].clone() * t + v_input
        };
        let c_fhe_key = {
            let t = if *input { &public_data.t_1.1 } else { &public_data.t_0.1 };
            cs_fhe_key[idx].clone() * t + v_fhe_key
        };
        ps.push(p);
        cs_input.push(c_input);
        cs_fhe_key.push(c_fhe_key);
    }
    // [TODO] Operation using the final preimage.
    vec![]
}
