use circuit::BenchCircuit;
use clap::{Parser, Subcommand};
use config::{RunBenchConfig, SimBenchNormConfig};
#[cfg(feature = "disk")]
use diamond_io::utils::calculate_tmp_size;
use diamond_io::{
    io::{
        eval::evaluate, obf::obfuscate, params::ObfuscationParams,
        utils::build_final_digits_circuit,
    },
    poly::{
        Poly, PolyElem, PolyParams,
        dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
            DCRTPolyUniformSampler, FinRingElem,
        },
        sampler::{DistType, PolyUniformSampler},
    },
    utils::{calculate_directory_size, init_tracing},
};
use num_traits::identities::One;

#[cfg(feature = "disk")]
use diamond_io::utils::log_mem;
use keccak_asm::Keccak256;
use num_bigint::BigUint;
use std::{
    fs::{self},
    path::{Path, PathBuf},
    sync::Arc,
};
use tracing::info;

pub mod circuit;
pub mod config;

/// Simple program to obfuscate and evaluate
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    RunBench {
        #[arg(short, long)]
        config: PathBuf,

        #[arg(short, long)]
        obf_dir: PathBuf,

        #[arg(short, long, default_value = "true")]
        verify: bool,

        #[arg(long)]
        add_num: usize,

        #[arg(long)]
        mul_num: usize,
    },
    SimBenchNorm {
        #[arg(short, long)]
        config: PathBuf,

        #[arg(short, long)]
        out_path: PathBuf,

        #[arg(long)]
        add_num: usize,

        #[arg(long)]
        mul_num: usize,
    },
}

#[tokio::main]
async fn main() {
    init_tracing();
    let command = Args::parse().command;
    match command {
        Commands::RunBench { config, obf_dir, verify, add_num, mul_num } => {
            let contents = fs::read_to_string(&config).unwrap();
            let dio_config: RunBenchConfig = toml::from_str(&contents).unwrap();
            let dir = Path::new(&obf_dir);
            if !dir.exists() {
                fs::create_dir(dir).unwrap();
            } else {
                // Clean it first to ensure no old files interfere
                fs::remove_dir_all(dir).unwrap();
                fs::create_dir(dir).unwrap();
            }
            let start_time = std::time::Instant::now();
            let params = DCRTPolyParams::new(
                dio_config.ring_dimension,
                dio_config.crt_depth,
                dio_config.crt_bits,
                dio_config.base_bits,
            );
            let log_base_q = params.modulus_digits();
            let switched_modulus = Arc::new(dio_config.switched_modulus);
            let public_circuit =
                BenchCircuit::new_add_mul(add_num, mul_num, log_base_q).as_poly_circuit();
            info!("Public circuit: {:?}", public_circuit.count_gates_by_type_vec());

            let obf_params = ObfuscationParams {
                params: params.clone(),
                switched_modulus,
                input_size: dio_config.input_size,
                level_width: dio_config.level_width,
                public_circuit,
                d: dio_config.d,
                encoding_sigma: dio_config.encoding_sigma,
                hardcoded_key_sigma: dio_config.hardcoded_key_sigma,
                p_sigma: dio_config.p_sigma,
                trapdoor_sigma: dio_config.trapdoor_sigma.unwrap_or_default(),
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
            >(obf_params.clone(), hardcoded_key.clone(), &mut rng, &obf_dir)
            .await;
            let obfuscation_time = start_time.elapsed();
            info!("Time to obfuscate: {:?}", obfuscation_time);

            let obf_size = calculate_directory_size(&obf_dir);
            info!("Obfuscation size: {obf_size} bytes");

            let input = dio_config.input;
            assert_eq!(input.len(), dio_config.input_size);

            let start_time = std::time::Instant::now();
            let output = evaluate::<
                DCRTPolyMatrix,
                DCRTPolyHashSampler<Keccak256>,
                DCRTPolyTrapdoorSampler,
                _,
            >(obf_params, &input, &obf_dir);
            let eval_time = start_time.elapsed();
            let total_time = obfuscation_time + eval_time;
            info!("Time for evaluation: {:?}", eval_time);
            info!("Total time: {:?}", total_time);
            if verify {
                let verify_circuit =
                    BenchCircuit::new_add_mul_verify(add_num, mul_num).as_poly_circuit();
                let input_coeffs: Vec<_> = input
                    .iter()
                    .cloned()
                    .map(|i| FinRingElem::constant(&params.modulus(), i as u64))
                    .collect();
                let input_poly = DCRTPoly::from_coeffs(&params, &input_coeffs);
                let eval = verify_circuit.eval(
                    &params,
                    &DCRTPoly::const_one(&params),
                    &[hardcoded_key, input_poly],
                );
                assert_eq!(eval.len(), 1);
                /*
                    Since we are computing b' - a' * t in the decryption part of the final circuit,
                    where a' = acc * a and b' = acc * b are the outputs of the public circuit,
                    b' - a' * t = acc * b - acc * a * t = acc * (a * t + e + [q/2] x - a*t) = acc * (e + [q/2] x) should hold,
                    where e is the LWE error and x is the hardcoded key.
                    If e = 0, it follows that b' - a' * t = acc * [q/2] x.
                */
                let half_q = FinRingElem::half_q(&params.modulus());
                for e in eval {
                    let expected_output = (DCRTPoly::from_const(&params, &half_q) * e)
                        .extract_bits_with_threshold(&params);
                    assert_eq!(output, expected_output);
                }
            }
        }
        Commands::SimBenchNorm { config, out_path, add_num, mul_num } => {
            let dio_config: SimBenchNormConfig =
                serde_json::from_reader(fs::File::open(&config).unwrap()).unwrap();
            let log_n = dio_config.log_ring_dim;
            let n = 2u32.pow(log_n);
            let max_crt_depth = dio_config.max_crt_depth;
            let crt_bits = dio_config.crt_bits;
            let base_bits = dio_config.base_bits;
            let params = DCRTPolyParams::new(n, max_crt_depth, crt_bits, base_bits);
            let log_q = params.modulus_bits();
            debug_assert_eq!(crt_bits * max_crt_depth, log_q);
            let log_base_q = params.modulus_digits();

            let public_circuit =
                BenchCircuit::new_add_mul(add_num, mul_num, log_base_q).as_poly_circuit();
            let a_rlwe_bar = DCRTPoly::const_max(&params);
            let b = DCRTPoly::const_max(&params);

            let a_decomposed_polys = a_rlwe_bar.decompose_base(&params);
            let b_decomposed_polys = b.decompose_base(&params);
            let final_circuit = build_final_digits_circuit::<DCRTPoly, DCRTPoly>(
                &a_decomposed_polys,
                &b_decomposed_polys,
                public_circuit,
            );

            let packed_input_norms = vec![BigUint::one(), params.modulus().as_ref().clone()];
            let norms = final_circuit.simulate_bgg_norm(
                params.ring_dimension(),
                params.base_bits(),
                packed_input_norms,
            );
            let norm_json = serde_json::to_string(&norms).unwrap();
            fs::write(out_path, norm_json.as_bytes()).unwrap()
        }
    }
}
