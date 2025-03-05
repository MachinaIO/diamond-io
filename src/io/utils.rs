use super::ObfuscationParams;
use crate::bgg::{sampler::*, BggPublicKey};
use crate::poly::{matrix::*, sampler::*, Poly, PolyParams};
use itertools::Itertools;
use std::marker::PhantomData;
use std::sync::Arc;

const TAG_R_0: &[u8] = b"R_0";
const TAG_R_1: &[u8] = b"R_1";
const TAG_A_FHE_BAR: &[u8] = b"A_FHE_BAR";
const TAG_BGG_PUBKEY_INPUT_PREFIX: &[u8] = b"BGG_PUBKEY_INPUT:";
const TAG_BGG_PUBKEY_FHEKEY_PREFIX: &[u8] = b"BGG_PUBKEY_FHEKY:";
const TAG_A_PRF: &[u8] = b"A_PRF:";

#[derive(Debug, Clone)]
pub struct PublicSampledData<S: PolyHashSampler<[u8; 32]>> {
    pub r_0: S::M,
    pub r_1: S::M,
    pub a_fhe_bar: S::M,
    pub pubkeys_input: Vec<Vec<BggPublicKey<S::M>>>,
    pub pubkeys_fhe_key: Vec<Vec<BggPublicKey<S::M>>>,
    pub t_0: (S::M, S::M),
    pub t_1: (S::M, S::M),
    pub a_prf: S::M,
    _s: PhantomData<S>,
}

impl<S: PolyHashSampler<[u8; 32]>> PublicSampledData<S> {
    pub fn sample(
        obf_params: &ObfuscationParams<S::M>,
        bgg_pubkey_sampler: &BGGPublicKeySampler<[u8; 32], S>,
        packed_input_size: usize,
        packed_output_size: usize,
    ) -> Self {
        let hash_sampler = &bgg_pubkey_sampler.sampler;
        let params = &obf_params.params;
        let r_0_bar = hash_sampler.sample_hash(params, TAG_R_0, 1, 1, DistType::BitDist);
        let r_1_bar = hash_sampler.sample_hash(params, TAG_R_1, 1, 1, DistType::BitDist);
        let one = S::M::identity(params, 1, None);
        let r_0 = r_0_bar.concat_diag(&[one.clone()]);
        let r_1 = r_1_bar.concat_diag(&[one.clone()]);
        let log_q = params.modulus_bits();
        let a_fhe_bar =
            hash_sampler.sample_hash(params, TAG_A_FHE_BAR, 2, 2 * log_q, DistType::FinRingDist);
        let pubkeys_input = (0..obf_params.input_size + 1)
            .map(|idx| {
                bgg_pubkey_sampler.sample(
                    params,
                    &[TAG_BGG_PUBKEY_INPUT_PREFIX, &idx.to_le_bytes()].concat(),
                    packed_input_size + 1,
                )
            })
            .collect_vec();
        let pubkeys_fhe_key = (0..obf_params.input_size + 1)
            .map(|idx| {
                bgg_pubkey_sampler.sample(
                    params,
                    &[TAG_BGG_PUBKEY_FHEKEY_PREFIX, &idx.to_le_bytes()].concat(),
                    2,
                )
            })
            .collect_vec();
        let identity_input = S::M::identity(params, packed_input_size + 1, None);
        let gadget_2 = S::M::gadget_matrix(params, 2);
        let identity_2 = S::M::identity(params, 2, None);
        let mut ts = vec![];
        for bit in 0..1 {
            let r = if bit == 0 { r_0.clone() } else { r_1.clone() };
            let rg = r * &gadget_2;
            let rg_decomposed = rg.decompose();
            let t_input = identity_input.clone().tensor(&rg_decomposed);
            let t_fhe_key = identity_2.clone().tensor(&rg_decomposed);
            ts.push((t_input, t_fhe_key));
        }
        let a_prf_raw = hash_sampler.sample_hash(
            &params,
            TAG_A_PRF,
            2,
            packed_output_size,
            DistType::FinRingDist,
        );
        let a_prf = a_prf_raw.modulus_switch(&obf_params.modulus_switch_params);
        Self {
            r_0,
            r_1,
            a_fhe_bar,
            pubkeys_input,
            pubkeys_fhe_key,
            t_0: ts[0].clone(),
            t_1: ts[1].clone(),
            a_prf,
            _s: PhantomData,
        }
    }
}
