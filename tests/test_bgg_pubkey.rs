use std::collections::HashMap;

use diamond_io::{
    bgg::{
        circuit::PolyCircuit, lut::public_lut::PublicLut, public_key::BggPubKeyPltEvaluator,
        sampler::BGGPublicKeySampler, BggPublicKey, DigitsToInt,
    },
    io::utils::build_final_digits_circuit,
    poly::{
        dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
            DCRTPolyUniformSampler,
        },
        enc::rlwe_encrypt,
        sampler::{DistType, PolyHashSampler, PolyUniformSampler},
        Poly, PolyParams,
    },
    utils::{init_tracing, log_mem},
};
use keccak_asm::Keccak256;
use rand::Rng;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};
use tracing::info;

#[test]
#[ignore]
fn test_build_final_step_circuit() {
    init_tracing();
    let params = DCRTPolyParams::default();
    let log_base_q = params.modulus_digits();
    let mut public_circuit = PolyCircuit::new();

    // inputs: BaseDecompose(ct), eval_input
    // outputs: BaseDecompose(ct) AND eval_input
    {
        let inputs = public_circuit.input((2 * log_base_q) + 1);
        let mut outputs = vec![];
        let eval_input = inputs[2 * log_base_q];
        for ct_input in inputs[0..2 * log_base_q].iter() {
            let muled = public_circuit.and_gate(*ct_input, eval_input);
            outputs.push(muled);
        }
        public_circuit.output(outputs);
    }

    let mut rng = rand::rng();
    let hash_key = rng.random::<[u8; 32]>();
    let d = 1;
    let sampler_uniform = DCRTPolyUniformSampler::new();
    let hash_sampler = DCRTPolyHashSampler::<Keccak256>::new();
    let hardcoded_key = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
    log_mem("Sampled hardcoded_key_matrix");

    let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
    log_mem("Sampled t_bar_matrix");

    let a_rlwe_bar =
        hash_sampler.sample_hash(&params, hash_key, "TEST_RLWE_A", 1, 1, DistType::FinRingDist);
    let hardcoded_key_sigma = 40615715852990820734.97011;

    let b = rlwe_encrypt(
        &params,
        &sampler_uniform,
        &t_bar_matrix,
        &a_rlwe_bar,
        &hardcoded_key,
        hardcoded_key_sigma,
    );
    log_mem("Generated RLWE ciphertext {a, b}");

    let a_decomposed = a_rlwe_bar.entry(0, 0).decompose_base(&params);
    let b_decomposed = b.entry(0, 0).decompose_base(&params);
    log_mem("Decomposed RLWE ciphertext into {BaseDecompose(a), BaseDecompose(b)}");

    let final_circuit = build_final_digits_circuit::<DCRTPoly, BggPublicKey<DCRTPolyMatrix>>(
        &a_decomposed,
        &b_decomposed,
        public_circuit.clone(),
    );
    log_mem("Computed final_circuit");

    let input_size = 1usize;
    let dim = params.ring_dimension() as usize;
    let log_base_q = params.modulus_digits();
    let packed_input_size = input_size.div_ceil(dim) + 1;

    let bgg_pubkey_sampler =
        BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(hash_key, d);
    // consider inserting t on plaintext
    let reveal_plaintexts = [vec![true; packed_input_size], vec![false; 1]].concat();
    let pubkeys = bgg_pubkey_sampler.sample(&params, b"BGG_PUBKEY_INPUT:", &reveal_plaintexts);
    log_mem("Sampled pubkeys");

    let eval_outputs = final_circuit.eval(
        &params,
        &pubkeys[0],
        &pubkeys[1..],
        None::<
            BggPubKeyPltEvaluator<
                DCRTPolyMatrix,
                DCRTPolyHashSampler<Keccak256>,
                DCRTPolyUniformSampler,
                DCRTPolyTrapdoorSampler,
            >,
        >,
    );
    log_mem("Evaluated outputs");

    let output_ints = eval_outputs
        .par_chunks(log_base_q)
        .map(|bits| BggPublicKey::digits_to_int(bits, &params))
        .collect::<Vec<_>>();
    log_mem("Converted outputs to integers");
    let _ = output_ints[0].concat_matrix(&output_ints[1..]);
    log_mem("Concatenated outputs into matrix");
}

// #[test]
// #[ignore]
// fn test_plt_build_final_step_circuit() {
//     init_tracing();
//     let params = DCRTPolyParams::default();
//     let log_base_q = params.modulus_digits();
//     let mut public_circuit = PolyCircuit::new();

//     /* bit repr => const int repr */
//     let mut f = HashMap::new();
//     f.insert(0, (DCRTPoly::from_const_int_lsb(&params, 0), DCRTPoly::const_int(&params, 0)));
//     f.insert(1, (DCRTPoly::from_const_int_lsb(&params, 1), DCRTPoly::const_int(&params, 1)));
//     f.insert(2, (DCRTPoly::from_const_int_lsb(&params, 2), DCRTPoly::const_int(&params, 2)));
//     f.insert(3, (DCRTPoly::from_const_int_lsb(&params, 3), DCRTPoly::const_int(&params, 3)));
//     f.insert(4, (DCRTPoly::from_const_int_lsb(&params, 4), DCRTPoly::const_int(&params, 4)));
//     f.insert(5, (DCRTPoly::from_const_int_lsb(&params, 5), DCRTPoly::const_int(&params, 5)));
//     f.insert(6, (DCRTPoly::from_const_int_lsb(&params, 6), DCRTPoly::const_int(&params, 6)));
//     f.insert(7, (DCRTPoly::from_const_int_lsb(&params, 7), DCRTPoly::const_int(&params, 7)));

//     let d = 1;
//     let lut = PublicLut::<DCRTPoly>::new(f);

//     // inputs: BaseDecompose(ct), eval_input
//     // outputs: BaseDecompose(ct) AND eval_input
//     let inputs = public_circuit.input((2 * log_base_q) + 1);
//     let mut outputs = vec![];
//     let eval_input = inputs[2 * log_base_q];
//     let plt_id = public_circuit.register_public_lookup(lut.clone());
//     let plt_gate = public_circuit.public_lookup_gate(eval_input, plt_id);
//     for ct_input in inputs[0..2 * log_base_q].iter() {
//         let muled = public_circuit.public_lookup_gate(*ct_input, plt_id);
//         outputs.push(muled);
//     }
//     outputs.push(plt_gate);
//     public_circuit.output(outputs);
//     info!("plt len {}", public_circuit.lookups.len());

//     info!("po {}", public_circuit.num_output());

//     let mut rng = rand::rng();
//     let hash_key = rng.random::<[u8; 32]>();

//     let sampler_uniform = DCRTPolyUniformSampler::new();
//     let hash_sampler = DCRTPolyHashSampler::<Keccak256>::new();
//     let hardcoded_key = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
//     log_mem("Sampled hardcoded_key_matrix");

//     let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
//     log_mem("Sampled t_bar_matrix");

//     let a_rlwe_bar =
//         hash_sampler.sample_hash(&params, hash_key, "TEST_RLWE_A", 1, 1, DistType::FinRingDist);
//     let hardcoded_key_sigma = 40615715852990820734.97011;

//     let b = rlwe_encrypt(
//         &params,
//         &sampler_uniform,
//         &t_bar_matrix,
//         &a_rlwe_bar,
//         &hardcoded_key,
//         hardcoded_key_sigma,
//     );
//     log_mem("Generated RLWE ciphertext {a, b}");

//     let a_decomposed = a_rlwe_bar.entry(0, 0).decompose_base(&params);
//     let b_decomposed = b.entry(0, 0).decompose_base(&params);
//     log_mem("Decomposed RLWE ciphertext into {BaseDecompose(a), BaseDecompose(b)}");

//     let final_circuit = build_final_digits_circuit::<DCRTPoly, BggPublicKey<DCRTPolyMatrix>>(
//         &a_decomposed,
//         &b_decomposed,
//         public_circuit.clone(),
//     );
//     log_mem("Computed final_circuit");
//     info!("{}", final_circuit.num_output());
//     info!("plt len {}", final_circuit.lookups.len());

//     let input_size = 1usize;
//     let dim = params.ring_dimension() as usize;
//     let log_base_q = params.modulus_digits();
//     let packed_input_size = input_size.div_ceil(dim) + 1;

//     let bgg_pubkey_sampler =
//         BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(hash_key, d);
//     let reveal_plaintexts = vec![true; packed_input_size + 1];
//     let pubkeys = bgg_pubkey_sampler.sample(&params, b"BGG_PUBKEY_INPUT:", &reveal_plaintexts);
//     log_mem("Sampled pubkeys");

//     let eval_outputs = final_circuit.eval(&params, &pubkeys[0], &pubkeys[1..], None);
//     log_mem("Evaluated outputs");

//     let output_ints = eval_outputs
//         .par_chunks(log_base_q)
//         .map(|bits| BggPublicKey::digits_to_int(bits, &params))
//         .collect::<Vec<_>>();
//     log_mem("Converted outputs to integers");
//     let _ = output_ints[0].concat_matrix(&output_ints[1..]);
//     log_mem("Concatenated outputs into matrix");
// }
