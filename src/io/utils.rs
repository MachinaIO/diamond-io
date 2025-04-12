use crate::{
    bgg::{
        circuit::{build_composite_circuit_from_public_and_fhe_dec, Evaluable, PolyCircuit},
        sampler::*,
        BggPublicKey,
    },
    poly::{sampler::*, Poly, PolyElem, PolyMatrix, PolyParams},
};
use itertools::Itertools;
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
        let d = obf_params.d;

        let r_0_bar = hash_sampler.sample_hash(params, TAG_R_0, d, d, DistType::BitDist);
        let r_1_bar = hash_sampler.sample_hash(params, TAG_R_1, d, d, DistType::BitDist);
        let one = S::M::identity(params, 1, None);
        let r_0 = r_0_bar.concat_diag(&[&one]);
        let r_1 = r_1_bar.concat_diag(&[&one]);
        let log_q = params.modulus_bits();
        let dim = params.ring_dimension() as usize;
        // input bits, poly of the RLWE key
        let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
        let packed_output_size = obf_params.public_circuit.num_output() / (2 * log_q);
        let a_rlwe_bar =
            hash_sampler.sample_hash(params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);
        let gadget_d_plus_1 = S::M::gadget_matrix(params, d + 1);
        let rgs: [<S as PolyHashSampler<[u8; 32]>>::M; 2] =
            [(r_0.clone() * &gadget_d_plus_1), (r_1.clone() * &gadget_d_plus_1)];

        let a_prf_raw = hash_sampler.sample_hash(
            params,
            TAG_A_PRF,
            d + 1,
            packed_output_size,
            DistType::FinRingDist,
        );
        let a_prf = a_prf_raw.modulus_switch(&obf_params.switched_modulus);
        Self {
            r_0,
            r_1,
            a_rlwe_bar,
            rgs,
            a_prf,
            packed_input_size,
            packed_output_size,
            _s: PhantomData,
        }
    }
}

pub fn build_final_bits_circuit<P: Poly, E: Evaluable<P>>(
    a_decomposed_polys: &[P],
    b_decomposed_polys: &[P],
    public_circuit: PolyCircuit<P>,
) -> PolyCircuit<P> {
    let log_q = a_decomposed_polys.len();
    debug_assert_eq!(b_decomposed_polys.len(), log_q);
    let packed_eval_input_size = public_circuit.num_input() - (2 * log_q);

    // circuit outputs the cipertext ct=(a,b) as a_bit_0, b_bit_0, a_bit_1, b_bit_1, ...
    let mut ct_output_circuit = PolyCircuit::new();
    {
        let inputs = ct_output_circuit.input(packed_eval_input_size);
        let circuit_id = ct_output_circuit.register_sub_circuit(public_circuit);
        let mut public_circuit_inputs = vec![];
        for poly in a_decomposed_polys.iter() {
            let bits = poly.coeffs().iter().map(|elem| elem.to_bit()).collect_vec();
            public_circuit_inputs.push(ct_output_circuit.const_bit_poly(&bits));
        }
        for poly in b_decomposed_polys.iter() {
            let bits = poly.coeffs().iter().map(|elem| elem.to_bit()).collect_vec();
            public_circuit_inputs.push(ct_output_circuit.const_bit_poly(&bits));
        }
        public_circuit_inputs.extend(inputs);
        assert_eq!(public_circuit_inputs.len(), 2 * log_q + packed_eval_input_size);
        let pc_outputs = ct_output_circuit.call_sub_circuit(circuit_id, &public_circuit_inputs);
        let mut outputs = Vec::with_capacity(pc_outputs.len());
        // n is the number of ciphertexts
        let n = pc_outputs.len() / (2 * log_q);
        for ct_idx in 0..n {
            let ct_offset = ct_idx * 2 * log_q;
            for bit_idx in 0..log_q {
                let a_index = ct_offset + bit_idx;
                let b_index = ct_offset + log_q + bit_idx;
                let a_bit = pc_outputs[a_index];
                let b_bit = pc_outputs[b_index];
                outputs.push(a_bit);
                outputs.push(b_bit);
            }
        }
        ct_output_circuit.output(outputs);
    }

    // actual
    let mut circuit = PolyCircuit::new();
    {
        let inputs = circuit.input(packed_eval_input_size + 1); // + 1 is for -t_bar
        let sub_circuit =
            build_composite_circuit_from_public_and_fhe_dec::<P, E>(ct_output_circuit, log_q);
        let circuit_id = circuit.register_sub_circuit(sub_circuit);
        let outputs = circuit.call_sub_circuit(circuit_id, &inputs);
        circuit.output(outputs);
    }
    circuit
}

#[cfg(test)]
#[cfg(feature = "test")]
mod test {
    use super::*;
    use crate::{
        bgg::BitToInt,
        poly::{
            dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler},
            enc::rlwe_encrypt,
            sampler::DistType,
        },
    };

    #[test]
    fn test_build_final_step_circuit() {
        // 1. Set up parameters
        let params = DCRTPolyParams::default();
        let log_q = params.modulus_bits();
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sigma = 3.0;

        // 2. Create a simple public circuit that takes 2*log_q inputs and outputs them directly
        let mut public_circuit = PolyCircuit::new();
        {
            let inputs = public_circuit.input(2 * log_q + 1);
            public_circuit.output(inputs[0..2 * log_q].to_vec());
        }

        // 3. Generate a random hardcoded key
        let hardcoded_key = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

        // 4. Generate RLWE ciphertext for the hardcoded key
        let a_rlwe_bar = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
        let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

        let b = rlwe_encrypt(
            &params,
            &sampler_uniform,
            &t_bar_matrix,
            &a_rlwe_bar,
            &hardcoded_key,
            sigma,
        );

        // 5. Decompose the ciphertext
        let a_decomposed = a_rlwe_bar.entry(0, 0).decompose_bits(&params);
        let b_decomposed = b.entry(0, 0).decompose_bits(&params);

        // 6. Build the final circuit with DCRTPoly as the Evaluable type
        let final_circuit = build_final_bits_circuit::<DCRTPoly, DCRTPoly>(
            &a_decomposed,
            &b_decomposed,
            public_circuit,
        );

        // 7. Evaluate the circuit
        let one = DCRTPoly::const_one(&params);

        let mut inputs = vec![one.clone()];
        inputs.push(-(t_bar_matrix.entry(0, 0)).clone());

        let circuit_outputs = final_circuit.eval(&params, &one, &inputs);
        assert_eq!(circuit_outputs.len(), log_q);

        // 8. Extract the output bits
        let output_ints = circuit_outputs
            .chunks(log_q)
            .map(|bits| DCRTPoly::bits_to_int(bits, &params))
            .collect_vec();
        assert_eq!(output_ints.len(), 1);
        let output_bits = output_ints
            .iter()
            .flat_map(|output| output.extract_bits_with_threshold(&params))
            .collect::<Vec<_>>();

        // 9. Verify that the output matches the hardcoded key bits
        assert_eq!(output_bits.len(), params.ring_dimension() as usize);
        assert_eq!(output_bits, hardcoded_key.entry(0, 0).to_bool_vec());
    }

    // #[test]
    // fn test_simulate_norm_final_bits_circuit() {
    //     // 1. Set up parameters
    //     let params = DCRTPolyParams::new(8192, 12, 51, 16);
    //     let log_q = params.modulus_bits();

    //     // 2. Create a simple public circuit that takes log_q inputs and outputs them directly
    //     let mut public_circuit = PolyCircuit::new();
    //     {
    //         let inputs = public_circuit.input((2 * log_q) + 1);
    //         public_circuit.output(inputs[0..(2 * log_q)].to_vec());
    //     }

    //     let a_rlwe_bar = DCRTPoly::const_max(&params);
    //     let enc_hardcoded_key = DCRTPoly::const_max(&params);

    //     let a_decomposed_polys = a_rlwe_bar.decompose_bits(&params);
    //     let b_decomposed_polys = enc_hardcoded_key.decompose_bits(&params);
    //     let final_circuit = build_final_bits_circuit::<DCRTPoly, DCRTPoly>(
    //         &a_decomposed_polys,
    //         &b_decomposed_polys,
    //         public_circuit,
    //     );

    //     let norms = final_circuit.simulate_bgg_norm(
    //         params.ring_dimension(),
    //         params.base_bits(),
    //         1 + params.ring_dimension() as usize,
    //     );
    //     let norm_json = serde_json::to_string(&norms).unwrap();
    //     use std::{fs::File, io::Write};
    //     let mut file = File::create("final_bits_norm.json").unwrap();
    //     file.write_all(norm_json.as_bytes()).unwrap();
    // }
}
