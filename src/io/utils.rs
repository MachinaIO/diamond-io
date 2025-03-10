use super::ObfuscationParams;
use crate::bgg::circuit::{build_circuit_ip_to_int, PolyCircuit};
use crate::bgg::{sampler::*, BggPublicKey, Evaluable};
use crate::poly::{matrix::*, sampler::*, Poly, PolyParams};
use itertools::Itertools;
use std::marker::PhantomData;

const TAG_R_0: &[u8] = b"R_0";
const TAG_R_1: &[u8] = b"R_1";
const TAG_A_RLWE_BAR: &[u8] = b"A_RLWE_BAR";
const TAG_BGG_PUBKEY_INPUT_PREFIX: &[u8] = b"BGG_PUBKEY_INPUT:";
const TAG_BGG_PUBKEY_FHEKEY_PREFIX: &[u8] = b"BGG_PUBKEY_FHEKY:";
const TAG_A_PRF: &[u8] = b"A_PRF:";

#[derive(Debug, Clone)]
pub struct PublicSampledData<S: PolyHashSampler<[u8; 32]>> {
    pub r_0: S::M,
    pub r_1: S::M,
    pub a_rlwe_bar: S::M,
    pub pubkeys: Vec<Vec<BggPublicKey<S::M>>>,
    // pub pubkeys_fhe_key: Vec<Vec<BggPublicKey<S::M>>>,
    pub ts: [S::M; 2],
    pub a_prf: S::M,
    pub packed_input_size: usize,
    pub packed_output_size: usize,
    _s: PhantomData<S>,
}

impl<S: PolyHashSampler<[u8; 32]>> PublicSampledData<S> {
    pub fn sample(
        obf_params: &ObfuscationParams<S::M>,
        bgg_pubkey_sampler: &BGGPublicKeySampler<[u8; 32], S>,
    ) -> Self {
        let hash_sampler = &bgg_pubkey_sampler.sampler;
        let params = &obf_params.params;
        let r_0_bar = hash_sampler.sample_hash(params, TAG_R_0, 1, 1, DistType::BitDist);
        let r_1_bar = hash_sampler.sample_hash(params, TAG_R_1, 1, 1, DistType::BitDist);
        let one = S::M::identity(params, 1, None);
        let r_0 = r_0_bar.concat_diag(&[one.clone()]);
        let r_1 = r_1_bar.concat_diag(&[one.clone()]);
        let log_q = params.modulus_bits();
        let dim = params.ring_dimension() as usize;
        // (bits of encrypted hardcoded key, input bits, poly of the FHE key)
        let packed_input_size = log_q + obf_params.input_size.div_ceil(dim) + 1;
        let packed_output_size = obf_params.public_circuit.num_output() / log_q;
        let a_rlwe_bar =
            hash_sampler.sample_hash(params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);
        // let reveal_plaintexts_fhe_key = vec![true; 2];
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();
        let pubkeys = (0..obf_params.input_size + 1)
            .map(|idx| {
                bgg_pubkey_sampler.sample(
                    params,
                    &[TAG_BGG_PUBKEY_INPUT_PREFIX, &idx.to_le_bytes()].concat(),
                    &reveal_plaintexts,
                )
            })
            .collect_vec();
        // let pubkeys_fhe_key = (0..obf_params.input_size + 1)
        //     .map(|idx| {
        //         bgg_pubkey_sampler.sample(
        //             params,
        //             &[TAG_BGG_PUBKEY_FHEKEY_PREFIX, &idx.to_le_bytes()].concat(),
        //             2,
        //         )
        //     })
        //     .collect_vec();
        let identity_input = S::M::identity(params, 1 + packed_input_size, None);
        let gadget_2 = S::M::gadget_matrix(params, 2);
        // let identity_2 = S::M::identity(params, 2, None);
        let mut ts = vec![];
        for bit in 0..2 {
            let r = if bit == 0 { r_0.clone() } else { r_1.clone() };
            let rg = r * &gadget_2;
            let rg_decomposed = rg.decompose();
            let t = identity_input.clone().tensor(&rg_decomposed);
            // let t_fhe_key = identity_2.clone().tensor(&rg_decomposed);
            ts.push(t);
        }
        let ts = ts.try_into().unwrap();
        let a_prf_raw = hash_sampler.sample_hash(
            params,
            TAG_A_PRF,
            2,
            packed_output_size,
            DistType::FinRingDist,
        );
        let a_prf = a_prf_raw.modulus_switch(&obf_params.switched_modulus);
        Self {
            r_0,
            r_1,
            a_rlwe_bar,
            pubkeys,
            ts,
            a_prf,
            packed_input_size,
            packed_output_size,
            _s: PhantomData,
        }
    }
}

pub fn build_final_step_circuit<P: Poly, E: Evaluable<P>>(
    params: &P::Params,
    a_decomposed_polys: &[P],
    public_circuit: PolyCircuit<P>,
) -> PolyCircuit<P> {
    let log_q = params.modulus_bits();
    debug_assert_eq!(a_decomposed_polys.len(), log_q);
    let packed_public_input_size = public_circuit.num_input();

    let mut ct_output_circuit = PolyCircuit::<P>::new();
    {
        let circuit_id = ct_output_circuit.register_sub_circuit(public_circuit);
        let inputs = ct_output_circuit.input(packed_public_input_size);
        let pc_outputs = ct_output_circuit.call_sub_circuit(circuit_id, &inputs);
        let mut outputs = vec![];
        for (idx, b_bit) in pc_outputs.iter().enumerate() {
            outputs.push(ct_output_circuit.const_scalar(a_decomposed_polys[idx].clone()));
            outputs.push(*b_bit);
            // let ct_bit = circuit.and_gate(*b_bit, inputs[packed_public_input_size]);
            // ct_bits.push(ct_bit);
        }
        ct_output_circuit.output(outputs);
    }
    let mut circuit = PolyCircuit::<P>::new();
    {
        let mut inputs = circuit.input(packed_public_input_size + 1);
        debug_assert_eq!(inputs.len(), packed_public_input_size + 1);
        let minus_one = circuit.const_minus_one_gate();
        inputs.push(minus_one);
        let sub_circuit = build_circuit_ip_to_int::<P, E>(params, ct_output_circuit, 2, log_q);
        let circuit_id = circuit.register_sub_circuit(sub_circuit);
        // debug_assert_eq!(public_outputs.len(), log_q);
        // let mut ct_bits = vec![];
        // for (idx, b_bit) in public_outputs.iter().enumerate() {
        //     ct_bits.push(circuit.const_scalar(a_decomposed_polys[idx].clone()));
        //     ct_bits.push(*b_bit);
        //     // let ct_bit = circuit.and_gate(*b_bit, inputs[packed_public_input_size]);
        //     // ct_bits.push(ct_bit);
        // }
        // inputs.extend(ct_bits);
        let outputs = circuit.call_sub_circuit(circuit_id, &inputs);
        println!("outputs len {:?}", outputs.len());
        circuit.output(outputs);
    }
    circuit
}
