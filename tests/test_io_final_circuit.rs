use diamond_io::{io::utils::build_final_digits_circuit, utils::init_tracing};
use keccak_asm::Keccak256;
use mxx::{
    bgg::{digits_to_int::DigitsToInt, public_key::BggPublicKey, sampler::BGGPublicKeySampler},
    circuit::PolyCircuit,
    lookup::simple_eval::SimpleBggPubKeyEvaluator,
    matrix::dcrt_poly::DCRTPolyMatrix,
    poly::{
        dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
        Poly, PolyParams,
    },
    rlwe_enc::rlwe_encrypt,
    sampler::{
        hash::DCRTPolyHashSampler, trapdoor::DCRTPolyTrapdoorSampler,
        uniform::DCRTPolyUniformSampler, DistType, PolyHashSampler, PolyUniformSampler,
    },
    utils::log_mem,
};
use rand::Rng;
use rayon::{iter::ParallelIterator, slice::ParallelSlice};

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
            SimpleBggPubKeyEvaluator<
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
