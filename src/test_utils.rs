use crate::{
    bgg::{circuit::PolyCircuit, lut::public_lut::PublicLut},
    io::{eval::evaluate, obf::obfuscate, params::ObfuscationParams},
    poly::{
        dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
            DCRTPolyUniformSampler, FinRingElem,
        },
        sampler::{DistType, PolyUniformSampler},
        Poly, PolyElem, PolyParams,
    },
    utils::{calculate_directory_size, init_tracing, log_mem},
};
use itertools::Itertools;
use keccak_asm::Keccak256;
use num_bigint::BigUint;
use num_traits::Num;
use rand::Rng;
use std::{collections::HashMap, fs, path::Path, sync::Arc};
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
    p_sigma: f64,
    hardcoded_key_sigma: f64,
    dir_path: &str,
) {
    init_tracing();
    let dir = Path::new(&dir_path);
    if !dir.exists() {
        fs::create_dir(dir).unwrap();
    } else {
        // Clean it first to ensure no old files interfere
        fs::remove_dir_all(dir).unwrap();
        fs::create_dir(dir).unwrap();
    }

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
        p_sigma,
        hardcoded_key_sigma,
        trapdoor_sigma: SIGMA,
    };

    let sampler_uniform = DCRTPolyUniformSampler::new();
    let mut rng = rand::rng();
    let hardcoded_key = sampler_uniform.sample_poly(&params, &DistType::BitDist);

    obfuscate::<
        DCRTPolyMatrix,
        DCRTPolyUniformSampler,
        DCRTPolyHashSampler<Keccak256>,
        DCRTPolyTrapdoorSampler,
        _,
        _,
    >(obf_params.clone(), hardcoded_key.clone(), &mut rng, &dir_path)
    .await;
    let obfuscation_time = start_time.elapsed();
    info!("Time to obfuscate: {:?}", obfuscation_time);

    let obf_size = calculate_directory_size(dir_path);
    log_mem(format!("Obfuscation size: {obf_size} bytes"));

    let bool_in = rng.random::<bool>();
    let mut input = vec![bool_in];
    input.append(&mut vec![false; input_size - 1]);

    let start_time = std::time::Instant::now();
    let output =
        evaluate::<DCRTPolyMatrix, DCRTPolyHashSampler<Keccak256>, DCRTPolyTrapdoorSampler, _>(
            obf_params, &input, &dir_path,
        );
    let eval_time = start_time.elapsed();
    info!("Time for evaluation: {:?}", eval_time);
    info!("Total time: {:?}", obfuscation_time + eval_time);

    let input_poly =
        DCRTPoly::from_const(&params, &FinRingElem::constant(&params.modulus(), bool_in as u64));
    assert_eq!(output, (hardcoded_key * input_poly).to_bool_vec());
}

pub async fn test_io_plt(
    ring_dim: u32,
    crt_depth: usize,
    crt_bits: usize,
    base_bits: u32,
    switched_modulus_str: &str,
    d: usize,
    input_size: usize,
    level_width: usize,
    p_sigma: f64,
    hardcoded_key_sigma: f64,
    dir_path: &str,
) {
    init_tracing();
    let dir = Path::new(&dir_path);
    if !dir.exists() {
        fs::create_dir(dir).unwrap();
    } else {
        // Clean it first to ensure no old files interfere
        fs::remove_dir_all(dir).unwrap();
        fs::create_dir(dir).unwrap();
    }

    let start_time = std::time::Instant::now();
    let params = DCRTPolyParams::new(ring_dim, crt_depth, crt_bits, base_bits);
    let log_base_q = params.modulus_digits();
    let switched_modulus = Arc::new(BigUint::from_str_radix(switched_modulus_str, 10).unwrap());
    let mut public_circuit = PolyCircuit::new();

    /* bit repr => const int repr */
    let mut f = HashMap::new();
    f.insert(
        0,
        (DCRTPoly::from_const_int_lsb(&params, 0), DCRTPoly::from_const_int_lsb(&params, 3)),
    );
    f.insert(
        1,
        (DCRTPoly::from_const_int_lsb(&params, 1), DCRTPoly::from_const_int_lsb(&params, 4)),
    );
    f.insert(
        2,
        (DCRTPoly::from_const_int_lsb(&params, 2), DCRTPoly::from_const_int_lsb(&params, 2)),
    );
    f.insert(
        3,
        (DCRTPoly::from_const_int_lsb(&params, 3), DCRTPoly::from_const_int_lsb(&params, 6)),
    );
    f.insert(
        4,
        (DCRTPoly::from_const_int_lsb(&params, 4), DCRTPoly::from_const_int_lsb(&params, 1)),
    );
    f.insert(
        5,
        (DCRTPoly::from_const_int_lsb(&params, 5), DCRTPoly::from_const_int_lsb(&params, 0)),
    );
    f.insert(
        6,
        (DCRTPoly::from_const_int_lsb(&params, 6), DCRTPoly::from_const_int_lsb(&params, 7)),
    );
    f.insert(
        7,
        (DCRTPoly::from_const_int_lsb(&params, 7), DCRTPoly::from_const_int_lsb(&params, 5)),
    );

    let lut = PublicLut::<DCRTPolyMatrix>::new::<
        DCRTPolyUniformSampler,
        DCRTPolyHashSampler<Keccak256>,
    >(&params, d, f, rand::random());

    // inputs: BaseDecompose(ct), eval_input
    // outputs: (eval_input PLT) * BaseDecompose(ct)
    let inputs = public_circuit.input((2 * log_base_q) + 1);
    let mut outputs = vec![];
    let eval_input = inputs[2 * log_base_q];
    let plt_id = public_circuit.register_public_lookup(lut.clone());
    let plt_out = public_circuit.public_lookup_gate(eval_input, plt_id);
    for ct_input in inputs[0..2 * log_base_q].iter() {
        let muled = public_circuit.mul_gate(*ct_input, plt_out);
        outputs.push(muled);
    }
    public_circuit.output(outputs);

    let obf_params = ObfuscationParams {
        params: params.clone(),
        switched_modulus,
        input_size,
        level_width,
        public_circuit: public_circuit.clone(),
        d,
        p_sigma,
        hardcoded_key_sigma,
        trapdoor_sigma: SIGMA,
    };

    let sampler_uniform = DCRTPolyUniformSampler::new();
    let mut rng = rand::rng();
    let hardcoded_key = sampler_uniform.sample_poly(&params, &DistType::BitDist);

    obfuscate::<
        DCRTPolyMatrix,
        DCRTPolyUniformSampler,
        DCRTPolyHashSampler<Keccak256>,
        DCRTPolyTrapdoorSampler,
        _,
        _,
    >(obf_params.clone(), hardcoded_key.clone(), &mut rng, &dir_path)
    .await;
    let obfuscation_time = start_time.elapsed();
    info!("Time to obfuscate: {:?}", obfuscation_time);

    let obf_size = calculate_directory_size(dir_path);
    log_mem(format!("Obfuscation size: {obf_size} bytes"));

    // 0,1,1(lsb) -> 4
    let input = vec![false, true, true];

    let start_time = std::time::Instant::now();
    let output =
        evaluate::<DCRTPolyMatrix, DCRTPolyHashSampler<Keccak256>, DCRTPolyTrapdoorSampler, _>(
            obf_params, &input, &dir_path,
        );
    let eval_time = start_time.elapsed();
    info!("output: {:?}", output);
    info!("Time for evaluation: {:?}", eval_time);
    info!("Total time: {:?}", obfuscation_time + eval_time);
    let output_poly = DCRTPoly::from_coeffs(
        &params,
        &output
            .iter()
            .map(|b| FinRingElem::from_bytes(&params.modulus(), &[*b as u8]))
            .collect_vec(),
    );
    let scale = DCRTPoly::from_const_int_lsb(&params, 7);
    // Public lookup for 0,1,1(lsb) => 1,1,1(lsb)
    // (1,1,1) * hardcoded key = output
    assert_eq!(output_poly.coeffs(), (hardcoded_key * scale).coeffs());
}
