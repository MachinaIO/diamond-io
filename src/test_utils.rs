use crate::{
    bgg::circuit::PolyCircuit,
    io::{obf::obfuscate, params::ObfuscationParams, Obfuscation},
    poly::{
        dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
            DCRTPolyUniformSampler, FinRingElem,
        },
        sampler::DistType,
        Poly, PolyElem, PolyParams,
    },
    utils::init_tracing,
};
use keccak_asm::Keccak256;
use num_bigint::BigUint;
use num_traits::Num;
use rand::Rng;
use std::sync::Arc;
use tracing::info;

const SIGMA: f64 = 4.578;

pub async fn test_io_common(
    ring_dim: u32,
    crt_depth: usize,
    crt_bits: usize,
    base_bits: u32,
    switched_modulus_str: &str,
    d: usize,
    input_size: usize,
    level_width: usize,
    encoding_sigma: f64,
    hardcoded_key_sigma: f64,
    p_sigma: f64,
    dir_path: &str,
) {
    init_tracing();
    let start_time = std::time::Instant::now();
    let params = DCRTPolyParams::new(ring_dim, crt_depth, crt_bits, base_bits);
    let log_base_q = params.modulus_digits();
    let switched_modulus = Arc::new(BigUint::from_str_radix(switched_modulus_str, 10).unwrap());
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

    let obf_params = ObfuscationParams {
        params: params.clone(),
        switched_modulus,
        input_size,
        level_width,
        public_circuit: public_circuit.clone(),
        d,
        encoding_sigma,
        hardcoded_key_sigma,
        p_sigma,
    };

    let sampler_uniform = DCRTPolyUniformSampler::new();
    let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
    let sampler_trapdoor = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
    let mut rng = rand::rng();
    let hardcoded_key = sampler_uniform.sample_poly(&params, &DistType::BitDist);

    obfuscate::<DCRTPolyMatrix, _, _, _, _, _>(
        obf_params.clone(),
        sampler_uniform,
        sampler_hash,
        sampler_trapdoor,
        hardcoded_key.clone(),
        &mut rng,
        &dir_path,
    )
    .await;
    let obfuscation_time = start_time.elapsed();
    info!("Time to obfuscate: {:?}", obfuscation_time);

    let bool_in = rng.random::<bool>();
    let mut input = vec![bool_in];
    input.append(&mut vec![false; input_size - 1]);
    let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
    let start_time = std::time::Instant::now();
    let obfuscation = Obfuscation::read_dir(&obf_params, dir_path);
    let load_time = start_time.elapsed();
    info!("Time to load obfuscation: {:?}", load_time);
    let start_time = std::time::Instant::now();
    let output = obfuscation.eval::<_, DCRTPolyTrapdoorSampler>(&obf_params, sampler_hash, &input);
    let eval_time = start_time.elapsed();
    info!("Time for evaluation: {:?}", eval_time);
    info!("Total time: {:?}", obfuscation_time + load_time + eval_time);
    let input_poly =
        DCRTPoly::from_const(&params, &FinRingElem::constant(&params.modulus(), bool_in as u64));
    assert_eq!(output, (hardcoded_key.clone() * input_poly.clone()).to_bool_vec());
}
