use circuit::BenchCircuit;
use clap::{Parser, Subcommand};
use config::Config;
#[cfg(feature = "disk")]
use diamond_io::utils::calculate_tmp_size;
use diamond_io::{
    io::{Obfuscation, obf::obfuscate, params::ObfuscationParams},
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
use keccak_asm::Keccak256;
use std::{
    fs::{self},
    path::{Path, PathBuf},
    sync::Arc,
};
use tokio;
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
}

#[tokio::main]
async fn main() {
    init_tracing();
    let command = Args::parse().command;
    match command {
        Commands::RunBench { config, obf_dir, verify, add_num, mul_num } => {
            let contents = fs::read_to_string(&config).unwrap();
            let dio_config: Config = toml::from_str(&contents).unwrap();
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

            #[cfg(feature = "disk")]
            let tmp_start_size = calculate_tmp_size();
            let start_time = std::time::Instant::now();
            let obfuscation = Obfuscation::read_dir(&obf_params, &obf_dir);
            let load_time = start_time.elapsed();
            info!("Time to load obfuscation: {:?}", load_time);
            #[cfg(feature = "disk")]
            let tmp_end_size = calculate_tmp_size();
            #[cfg(feature = "disk")]
            log_mem(format!(
                "Obfuscation size after loading: {} bytes",
                tmp_end_size - tmp_start_size
            ));
            let start_time = std::time::Instant::now();
            let output = obfuscation
                .eval::<DCRTPolyHashSampler<Keccak256>, DCRTPolyTrapdoorSampler>(
                    obf_params, &input,
                );
            let eval_time = start_time.elapsed();
            let total_time = obfuscation_time + load_time + eval_time;
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
                /*
                    Since we are computing b' - a' * t in the decryption part of the final circuit,
                    where a' = acc * a and b' = acc * b are the outputs of the public circuit,
                    b' - a' * t = acc * b - acc * a * t = acc * (a * t + e + [q/2] x - a*t) = acc * (e + [q/2] x) should hold,
                    where e is the LWE error and x is the hardcoded key.
                    If e = 0, it follows that b' - a' * t = acc * [q/2] x.
                */
                let half_q = FinRingElem::half_q(&params.modulus());
                let multiplied_hardcoded_key =
                    hardcoded_key * DCRTPoly::from_const(&params, &half_q);
                let eval = verify_circuit.eval(
                    &params,
                    &DCRTPoly::const_one(&params),
                    &[multiplied_hardcoded_key, input_poly],
                );
                assert_eq!(eval.len(), 1);
                for e in eval {
                    let decompose_poly = e.extract_bits_with_threshold(&params);
                    assert_eq!(output, decompose_poly);
                }
            }
        }
    }
}
