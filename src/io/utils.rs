use crate::{
    bgg::{
        circuit::{build_circuit_ip_priv_and_pub_outputs, Evaluable, PolyCircuit},
        sampler::*,
        BggPublicKey,
    },
    poly::{matrix::*, sampler::*, Poly, PolyElem, PolyParams},
};
use itertools::Itertools;
use std::{marker::PhantomData, ops::Mul};
use tracing::info;

use super::params::ObfuscationParams;

const TAG_R_0: &[u8] = b"R_0";
const TAG_R_1: &[u8] = b"R_1";
const TAG_A_RLWE_BAR: &[u8] = b"A_RLWE_BAR";
const TAG_BGG_PUBKEY_INPUT_PREFIX: &[u8] = b"BGG_PUBKEY_INPUT:";
const _TAG_BGG_PUBKEY_FHEKEY_PREFIX: &[u8] = b"BGG_PUBKEY_FHEKY:";
const TAG_A_PRF: &[u8] = b"A_PRF:";

#[derive(Debug, Clone)]
pub struct PublicSampledData<S: PolyHashSampler<[u8; 32]>> {
    pub r_0: S::M,
    pub r_1: S::M,
    pub a_rlwe_bar: S::M,
    pub pubkeys: Vec<Vec<BggPublicKey<S::M>>>,
    // pub pubkeys_fhe_key: Vec<Vec<BggPublicKey<S::M>>>,
    pub rgs: [S::M; 2],
    pub a_prf: S::M,
    pub packed_input_size: usize,
    pub packed_output_size: usize,
    _s: PhantomData<S>,
}

impl<S: PolyHashSampler<[u8; 32]>> PublicSampledData<S>
where
    for<'a> &'a S::M: Mul<&'a S::M, Output = S::M>,
{
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
        // (bits of encrypted hardcoded key, input bits, poly of the FHE key)
        let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
        let packed_output_size = obf_params.public_circuit.num_output() / log_q;
        let a_rlwe_bar =
            hash_sampler.sample_hash(params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);
        // let reveal_plaintexts_fhe_key = vec![true; 2];
        #[cfg(feature = "test")]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
        #[cfg(not(feature = "test"))]
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
        // let identity_input = S::M::identity(params, 1 + packed_input_size, None);
        let gadget_d1 = S::M::gadget_matrix(params, d + 1);
        // let identity_2 = S::M::identity(params, 2, None);
        let rgs: [<S as PolyHashSampler<[u8; 32]>>::M; 2] =
            [(&r_0 * &gadget_d1), (&r_1 * &gadget_d1)];

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
            pubkeys,
            rgs,
            a_prf,
            packed_input_size,
            packed_output_size,
            _s: PhantomData,
        }
    }
}

pub fn build_final_bits_circuit<P: Poly, E: Evaluable>(
    a_decomposed_polys: &[P],
    b_decomposed_polys: &[P],
    public_circuit: PolyCircuit,
) -> PolyCircuit {
    let log_q = a_decomposed_polys.len();
    debug_assert_eq!(b_decomposed_polys.len(), log_q);
    let packed_eval_input_size = public_circuit.num_input() - log_q;

    // circuit outputs the cipertext
    let mut ct_output_circuit = PolyCircuit::new();
    {
        let inputs = ct_output_circuit.input(packed_eval_input_size);
        let circuit_id = ct_output_circuit.register_sub_circuit(public_circuit);
        let mut public_circuit_inputs = vec![];
        for poly in b_decomposed_polys.iter() {
            let bits = poly.coeffs().iter().map(|elem| elem.to_bit()).collect_vec();
            public_circuit_inputs.push(ct_output_circuit.const_bit_poly(&bits));
        }
        public_circuit_inputs.extend(inputs);
        let pc_outputs = ct_output_circuit.call_sub_circuit(circuit_id, &public_circuit_inputs);
        let mut outputs = Vec::with_capacity(pc_outputs.len() * 2);
        for (idx, b_bit) in pc_outputs.into_iter().enumerate() {
            let bits =
                a_decomposed_polys[idx].coeffs().iter().map(|elem| elem.to_bit()).collect_vec();
            outputs.push(ct_output_circuit.const_bit_poly(&bits));
            outputs.push(b_bit);
        }
        ct_output_circuit.output(outputs);
    }

    // actual
    let mut circuit = PolyCircuit::new();
    {
        let mut inputs = circuit.input(packed_eval_input_size + 1);
        let minus_one = circuit.const_minus_one_gate();
        inputs.push(minus_one);
        let sub_circuit = build_circuit_ip_priv_and_pub_outputs::<E>(ct_output_circuit, 2);
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
            dcrt::{
                DCRTPoly, DCRTPolyHashSampler, DCRTPolyParams, DCRTPolyUniformSampler, FinRingElem,
            },
            element::PolyElem,
            sampler::DistType,
        },
    };
    use keccak_asm::Keccak256;
    use num_bigint::BigUint;

    #[test]
    fn test_build_final_step_circuit() {
        // 1. Set up parameters
        let params = DCRTPolyParams::default();
        let log_q = params.modulus_bits();

        // 2. Create a simple public circuit that takes log_q inputs and outputs them directly
        let mut public_circuit = PolyCircuit::new();
        {
            let inputs = public_circuit.input(log_q + 1);
            public_circuit.output(inputs[0..log_q].to_vec());
        }

        // 3. Generate a random hardcoded key (similar to obf.rs lines 53-65)
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let hardcoded_key = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

        // 4. Get the public polynomial a from PublicSampledData's sample function
        let hash_sampler = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let a_rlwe_bar =
            hash_sampler.sample_hash(&params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);

        // 5. Generate RLWE ciphertext for the hardcoded key
        let t_bar = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let e = sampler_uniform.sample_uniform(&params, 1, 1, DistType::GaussDist { sigma: 0.0 });

        // Create a scale value (half of q)
        let modulus = params.modulus();
        let half_q = FinRingElem::half_q(&modulus.clone());
        let scale = DCRTPoly::from_const(&params, &half_q);
        let enc_hardcoded_key =
            t_bar.clone() * &a_rlwe_bar + &e - &(hardcoded_key.clone() * &scale);

        // 6. Decompose the ciphertext
        let enc_hardcoded_key_polys =
            enc_hardcoded_key.get_column_matrix_decompose(0).get_column(0);
        // println!("enc_hardcoded_key_polys {}", enc_hardcoded_key_polys);

        // 7. Build the final step circuit with DCRTPoly as the Evaluable type
        let a_decomposed_polys = a_rlwe_bar.get_column_matrix_decompose(0).get_column(0);
        let final_circuit = build_final_bits_circuit::<DCRTPoly, DCRTPoly>(
            &a_decomposed_polys,
            &enc_hardcoded_key_polys,
            public_circuit,
        );

        // 8. Evaluate the circuit with the decomposed ciphertext
        let one = DCRTPoly::const_one(&params);
        let mut inputs = vec![one.clone()];
        inputs.push(t_bar.entry(0, 0).clone());
        let circuit_outputs = final_circuit.eval(&params, one, &inputs);
        // 9. Extract the hardcoded key bits
        let hardcoded_key_bits = hardcoded_key
            .entry(0, 0)
            .coeffs()
            .iter()
            .map(|elem| elem.value() != &BigUint::from(0u8))
            .collect::<Vec<_>>();

        // 10. Extract the output bits
        let output_ints: Vec<DCRTPoly> = circuit_outputs
            .chunks(log_q)
            .map(|bits| DCRTPoly::bits_to_int(bits, &params))
            .collect_vec();
        let output_bits =
            output_ints.iter().flat_map(|output| output.extract_highest_bits()).collect::<Vec<_>>();
        // 11. Verify that the output matches the hardcoded key bits
        assert_eq!(output_bits.len(), hardcoded_key_bits.len());
        for (i, (output_bit, key_bit)) in
            output_bits.iter().zip(hardcoded_key_bits.iter()).enumerate()
        {
            assert_eq!(output_bit, key_bit, "Bit mismatch at position {}", i);
        }
    }

    #[test]
    fn test_simulate_norm_final_bits_circuit() {
        // 1. Set up parameters
        let params = DCRTPolyParams::new(4096, 12, 51);
        let log_q = params.modulus_bits();

        // 2. Create a simple public circuit that takes log_q inputs and outputs them directly
        let mut public_circuit = PolyCircuit::new();
        {
            let inputs = public_circuit.input(log_q + 1);
            public_circuit.output(inputs[0..log_q].to_vec());
        }

        let a_rlwe_bar = DCRTPoly::const_max(&params);
        let enc_hardcoded_key = DCRTPoly::const_max(&params);

        let a_decomposed_polys = a_rlwe_bar.decompose(&params);
        let b_decomposed_polys = enc_hardcoded_key.decompose(&params);
        let final_circuit = build_final_bits_circuit::<DCRTPoly, DCRTPoly>(
            &a_decomposed_polys,
            &b_decomposed_polys,
            public_circuit,
        );

        let norms = final_circuit
            .simulate_bgg_norm(params.ring_dimension(), params.ring_dimension() as usize + 1);
        let norm_json = serde_json::to_string(&norms).unwrap();
        println!("norms: {}", norm_json);
        use std::{fs::File, io::Write};
        let mut file = File::create("final_bits_norm.json").unwrap();
        file.write_all(norm_json.as_bytes()).unwrap();
    }
}
