use super::ObfuscationParams;
use crate::{
    bgg::{
        circuit::{build_circuit_ip_to_int, PolyCircuit},
        sampler::*,
        BggPublicKey, Evaluable,
    },
    parallel_iter,
    poly::{matrix::*, sampler::*, Poly, PolyParams},
};
use itertools::Itertools;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::{marker::PhantomData, ops::Mul};
use tracing::info;

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
    pub rgs_decomposed: [S::M; 2],
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

        let one = S::M::identity(params, 1, None);

        let log_q = params.modulus_bits();
        let dim = params.ring_dimension() as usize;
        // (bits of encrypted hardcoded key, input bits, poly of the FHE key)
        let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
        let packed_output_size = obf_params.public_circuit.num_output() / log_q;
        let a_rlwe_bar =
            hash_sampler.sample_hash(params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);
        // let reveal_plaintexts_fhe_key = vec![true; 2];
        #[cfg(test)]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![true; 1]].concat();
        #[cfg(not(test))]
        let reveal_plaintexts = [vec![true; packed_input_size - 1], vec![false; 1]].concat();
        let pubkeys = (0..obf_params.input_size + 1)
            .map(|idx| {
                info!("try pubkey 1 computed");
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
        info!("pubkeys computed");
        let gadget_2 = S::M::gadget_matrix(params, 2);
        info!("gadget_2 computed");
        let r_0 =
            hash_sampler.sample_hash(params, TAG_R_0, 1, 1, DistType::BitDist).concat_diag(&[&one]);
        let r_1 =
            hash_sampler.sample_hash(params, TAG_R_1, 1, 1, DistType::BitDist).concat_diag(&[&one]);
        let rgs_decomposed: [<S as PolyHashSampler<[u8; 32]>>::M; 2] = {
            parallel_iter!(0..2)
                .map(|i| {
                    if i == 0 {
                        (&r_0 * &gadget_2).decompose()
                    } else {
                        (&r_1 * &gadget_2).decompose()
                    }
                })
                .collect::<Vec<_>>()
                .try_into()
                .unwrap()
        };
        info!("rgs_decomposed computed");
        let a_prf_raw = hash_sampler.sample_hash(
            params,
            TAG_A_PRF,
            2,
            packed_output_size,
            DistType::FinRingDist,
        );
        info!("a_prf_raw computed");
        let a_prf = a_prf_raw.modulus_switch(&obf_params.switched_modulus);
        info!("modulus_switch");
        Self {
            r_0,
            r_1,
            a_rlwe_bar,
            pubkeys,
            rgs_decomposed,
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
    b_decomposed_polys: &[P],
    public_circuit: PolyCircuit<P>,
) -> PolyCircuit<P> {
    let log_q = params.modulus_bits();
    debug_assert_eq!(a_decomposed_polys.len(), log_q);
    debug_assert_eq!(b_decomposed_polys.len(), log_q);
    let packed_eval_input_size = public_circuit.num_input() - log_q;

    // circuit outputs the ciper text
    let mut ct_output_circuit = PolyCircuit::<P>::new();
    {
        let inputs = ct_output_circuit.input(packed_eval_input_size);
        let circuit_id = ct_output_circuit.register_sub_circuit(public_circuit);
        let mut public_circuit_inputs = vec![];
        for poly in b_decomposed_polys.iter() {
            public_circuit_inputs.push(ct_output_circuit.const_scalar(poly.clone()));
        }
        public_circuit_inputs.extend(inputs);
        let pc_outputs = ct_output_circuit.call_sub_circuit(circuit_id, &public_circuit_inputs);
        let mut outputs = Vec::with_capacity(pc_outputs.len() * 2);
        for (idx, b_bit) in pc_outputs.into_iter().enumerate() {
            outputs.push(ct_output_circuit.const_scalar(a_decomposed_polys[idx].clone()));
            outputs.push(b_bit);
            // let ct_bit = circuit.and_gate(*b_bit, inputs[packed_public_input_size]);
            // ct_bits.push(ct_bit);
        }
        ct_output_circuit.output(outputs);
    }

    let mut circuit = PolyCircuit::<P>::new();
    {
        let mut inputs = circuit.input(packed_eval_input_size + 1);
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
        circuit.output(outputs);
    }
    circuit
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        bgg::{
            sampler::{BGGEncodingSampler, BGGPublicKeySampler},
            BggEncoding,
        },
        poly::{
            dcrt::{
                DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams,
                DCRTPolyUniformSampler, FinRingElem,
            },
            element::PolyElem,
            sampler::DistType,
        },
        utils::create_bit_random_poly,
    };
    use keccak_asm::Keccak256;
    use num_bigint::BigUint;
    use std::sync::Arc;

    #[test]
    fn test_build_final_step_circuit() {
        // 1. Set up parameters
        let params = DCRTPolyParams::default();
        let log_q = params.modulus_bits();

        // 2. Create a simple public circuit that takes log_q inputs and outputs them directly
        let mut public_circuit = PolyCircuit::<DCRTPoly>::new();
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
        let enc_hardcoded_key_polys = enc_hardcoded_key.decompose().get_column(0);

        // 7. Build the final step circuit with DCRTPoly as the Evaluable type
        let a_decomposed_polys = a_rlwe_bar.decompose().get_column(0);
        let final_circuit = build_final_step_circuit::<DCRTPoly, DCRTPoly>(
            &params,
            &a_decomposed_polys,
            &enc_hardcoded_key_polys,
            public_circuit,
        );

        // 8. Evaluate the circuit with the decomposed ciphertext
        let one = DCRTPoly::const_one(&params);
        let mut inputs = vec![one.clone()];
        inputs.push(t_bar.entry(0, 0).clone());
        let outputs = final_circuit.eval(&(), one, &inputs);

        // 9. Extract the hardcoded key bits
        let hardcoded_key_bits = hardcoded_key
            .entry(0, 0)
            .coeffs()
            .iter()
            .map(|elem| elem.value() != &BigUint::from(0u8))
            .collect::<Vec<_>>();

        // 10. Extract the output bits
        let output_bits =
            outputs.iter().flat_map(|output| output.extract_highest_bits()).collect::<Vec<_>>();
        // 11. Verify that the output matches the hardcoded key bits
        assert_eq!(output_bits.len(), hardcoded_key_bits.len());
        for (i, (output_bit, key_bit)) in
            output_bits.iter().zip(hardcoded_key_bits.iter()).enumerate()
        {
            assert_eq!(output_bit, key_bit, "Bit mismatch at position {}", i);
        }
    }

    #[test]
    fn test_build_final_step_circuit_with_bgg_encoding() {
        // 1. Set up parameters
        let params = DCRTPolyParams::default();
        let log_q = params.modulus_bits();

        // 2. Create a simple public circuit that takes log_q inputs and outputs them directly
        let mut public_circuit = PolyCircuit::<DCRTPoly>::new();
        {
            let inputs = public_circuit.input(log_q + 1);
            public_circuit.output(inputs[0..log_q].to_vec());
        }

        // 3. Generate a random hardcoded key (similar to obf.rs lines 53-65)
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let hardcoded_key = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

        // 4. Set up samplers for BggEncoding
        let hash_key = [0u8; 32];
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(hash_key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler.clone());
        let uniform_sampler = Arc::new(sampler_uniform);

        // 5. Get the public polynomial a from PublicSampledData's sample function
        let a_rlwe_bar =
            hash_sampler.sample_hash(&params, TAG_A_RLWE_BAR, 1, 1, DistType::FinRingDist);

        // 6. Generate RLWE ciphertext for the hardcoded key
        let t_bar = uniform_sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let e = uniform_sampler.sample_uniform(&params, 1, 1, DistType::GaussDist { sigma: 0.0 });

        // Create a scale value (half of q)
        let modulus = params.modulus();
        let half_q = FinRingElem::half_q(&modulus.clone());
        let scale = DCRTPoly::from_const(&params, &half_q);
        let enc_hardcoded_key =
            t_bar.clone() * &a_rlwe_bar + &e - &(hardcoded_key.clone() * &scale);

        // 7. Decompose the ciphertext
        let enc_hardcoded_key_polys = enc_hardcoded_key.decompose().get_column(0);

        // 8. Build the final step circuit with BggEncoding as the Evaluable type
        let a_decomposed_polys = a_rlwe_bar.decompose().get_column(0);
        let final_circuit_pubkey = build_final_step_circuit::<DCRTPoly, BggPublicKey<DCRTPolyMatrix>>(
            &params,
            &a_decomposed_polys,
            &enc_hardcoded_key_polys,
            public_circuit.clone(),
        );
        let final_circuit_encoding =
            build_final_step_circuit::<DCRTPoly, BggEncoding<DCRTPolyMatrix>>(
                &params,
                &a_decomposed_polys,
                &enc_hardcoded_key_polys,
                public_circuit,
            );

        // 9. Create BggEncoding instances for the inputs
        // First, create a secret key for the BGGEncodingSampler
        let secret = create_bit_random_poly(&params);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create public keys for the encodings
        let reveal_plaintexts = vec![true; 1 + 1]; // +1 for t_bar
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create plaintexts from the decomposed ciphertext
        let mut plaintexts = vec![DCRTPoly::const_one(&params)];
        plaintexts.push(t_bar.entry(0, 0).clone());

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);

        // 10. Evaluate the circuit with the BggPublicKey/BggEncoding inputs
        let outputs_pubkeys = final_circuit_pubkey.eval(&params, pubkeys[0].clone(), &pubkeys[1..]);
        let outputs_encodings =
            final_circuit_encoding.eval(&params, encodings[0].clone(), &encodings[1..]);
        assert_eq!(
            outputs_pubkeys[0].concat_matrix(&outputs_pubkeys[1..]),
            outputs_encodings[0].pubkey.concat_matrix(
                &outputs_encodings[1..].iter().map(|enc| enc.pubkey.clone()).collect_vec()
            )
        );

        // 11. Extract the hardcoded key bits
        let hardcoded_key_bits = hardcoded_key
            .entry(0, 0)
            .coeffs()
            .iter()
            .map(|elem| elem.value() != &BigUint::from(0u8))
            .collect::<Vec<_>>();

        // 12. Extract the output bits from the BggEncoding outputs
        let output_bits = outputs_encodings
            .iter()
            .flat_map(|output| {
                // Extract the plaintext from the BggEncoding and get its highest bits
                output.plaintext.as_ref().unwrap().extract_highest_bits()
            })
            .collect::<Vec<_>>();

        // 13. Verify that the output matches the hardcoded key bits
        assert_eq!(output_bits.len(), hardcoded_key_bits.len());
        for (i, (output_bit, key_bit)) in
            output_bits.iter().zip(hardcoded_key_bits.iter()).enumerate()
        {
            assert_eq!(output_bit, key_bit, "Bit mismatch at position {}", i);
        }

        let expected_vector = {
            let plaintext = DCRTPolyMatrix::from_poly_vec_row(
                &params,
                vec![outputs_encodings[0].plaintext.clone().unwrap()],
            );
            bgg_encoding_sampler.secret_vec *
                (outputs_encodings[0].pubkey.matrix.clone() -
                    plaintext.tensor(&DCRTPolyMatrix::gadget_matrix(&params, 2)))
        };
        assert_eq!(outputs_encodings[0].vector, expected_vector);
    }
}
