use crate::{
    bgg::{sampler::*, BggPublicKey},
    poly::{matrix::*, sampler::*, Poly, PolyParams},
};
use std::marker::PhantomData;

use super::params::ObfuscationParams;

const TAG_R_0: &[u8] = b"R_0";
const TAG_R_1: &[u8] = b"R_1";
const TAG_A_RLWE_BAR: &[u8] = b"A_RLWE_BAR";
const _TAG_BGG_PUBKEY_FHEKEY_PREFIX: &[u8] = b"BGG_PUBKEY_FHEKY:";
const TAG_A_PRF: &[u8] = b"A_PRF:";
pub const TAG_BGG_PUBKEY_INPUT_PREFIX: &[u8] = b"BGG_PUBKEY_INPUT:";

pub fn sample_public_key_by_idx<K: AsRef<[u8]>, S>(
    sampler: &BGGPublicKeySampler<K, S>,
    params: &<<<S as PolyHashSampler<K>>::M as PolyMatrix>::P as Poly>::Params,
    idx: usize,
    reveal_plaintexts: &[bool],
) -> Vec<BggPublicKey<<S as PolyHashSampler<K>>::M>>
where
    S: PolyHashSampler<K>,
{
    sampler.sample(
        params,
        &[TAG_BGG_PUBKEY_INPUT_PREFIX, &(idx as u64).to_le_bytes()].concat(),
        reveal_plaintexts,
    )
}

#[derive(Debug, Clone)]
pub struct PublicSampledData<S: PolyHashSampler<[u8; 32]>> {
    pub r_0: S::M,
    pub r_1: S::M,
    pub a_rlwe_bar: S::M,
    pub rgs: [S::M; 2],
    pub a_prf: S::M,
    pub packed_input_size: usize,
    _s: PhantomData<S>,
}

impl<S: PolyHashSampler<[u8; 32]>> PublicSampledData<S> {
    pub fn sample(
        obf_params: &ObfuscationParams<S::M>,
        bgg_pubkey_sampler: &BGGPublicKeySampler<[u8; 32], S>,
    ) -> Self {
        let hash_sampler = &bgg_pubkey_sampler.sampler;
        let params = &obf_params.params;
        let d = obf_params.d;

        let r_0_bar = hash_sampler.sample_hash(params, TAG_R_0, d, d, DistType::BitDist);
        let r_1_bar = hash_sampler.sample_hash(params, TAG_R_1, d, d, DistType::BitDist);
        let one = S::M::identity(params, 1, None);
        let r_0 = r_0_bar.concat_diag(&[&one]);
        let r_1 = r_1_bar.concat_diag(&[&one]);
        let dim = params.ring_dimension() as usize;
        let packed_input_size = obf_params.input_size.div_ceil(dim) + 1; // input bits + poly of the RLWE key
        let a_rlwe_bar =
            hash_sampler.sample_hash(params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);
        let gadget_d_plus_1 = S::M::gadget_matrix(params, d + 1);
        let rgs: [<S as PolyHashSampler<[u8; 32]>>::M; 2] =
            [(r_0.clone() * &gadget_d_plus_1), (r_1.clone() * &gadget_d_plus_1)];

        let a_prf_raw =
            hash_sampler.sample_hash(params, TAG_A_PRF, d + 1, 1, DistType::FinRingDist);
        let a_prf = a_prf_raw.modulus_switch(&obf_params.switched_modulus);
        Self { r_0, r_1, a_rlwe_bar, rgs, a_prf, packed_input_size, _s: PhantomData }
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod test {
    // #[test]
    // fn test_simulate_norm_final_bits_circuit() {
    //     // 1. Set up parameters
    //     let params = DCRTPolyParams::new(4096, 12, 51);
    //     let log_q = params.modulus_bits();

    //     // 2. Create a simple public circuit that takes log_q inputs and outputs them directly
    //     let mut public_circuit = PolyCircuit::new();
    //     {
    //         let inputs = public_circuit.input(log_q + 1);
    //         public_circuit.output(inputs[0..log_q].to_vec());
    //     }

    //     let a_rlwe_bar = DCRTPoly::const_max(&params);
    //     let enc_hardcoded_key = DCRTPoly::const_max(&params);

    //     let a_decomposed_polys = a_rlwe_bar.decompose(&params);
    //     let b_decomposed_polys = enc_hardcoded_key.decompose(&params);
    //     let final_circuit = build_final_bits_circuit::<DCRTPoly, DCRTPoly>(
    //         &a_decomposed_polys,
    //         &b_decomposed_polys,
    //         public_circuit,
    //     );

    //     let norms = final_circuit
    //         .simulate_bgg_norm(params.ring_dimension(), params.ring_dimension() as usize + 1);
    //     let norm_json = serde_json::to_string(&norms).unwrap();
    //     println!("norms: {}", norm_json);
    //     use std::{fs::File, io::Write};
    //     let mut file = File::create("final_bits_norm.json").unwrap();
    //     file.write_all(norm_json.as_bytes()).unwrap();
    // }
}
