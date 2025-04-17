use circuit::BenchCircuit;
use clap::{Parser, Subcommand};
use config::Config;
use diamond_io::{
    io::{obf::obfuscate, params::ObfuscationParams},
    poly::{
        Poly, PolyElem, PolyParams,
        dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
            DCRTPolyUniformSampler, FinRingElem,
        },
        sampler::{DistType, PolyUniformSampler},
    },
    utils::init_tracing,
};
use keccak_asm::Keccak256;
use std::{
    fs::{self},
    path::PathBuf,
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
    Run {
        #[arg(short, long)]
        config: PathBuf,

        #[arg(short, long, default_value = "true")]
        verify: bool,

        #[arg(long)]
        add_num: usize,

        #[arg(long)]
        mul_num: usize,
    },
}

fn main() {
    init_tracing();
    let command = Args::parse().command;
    match command {
        Commands::Run { config, verify, add_num, mul_num } => {
            let contents = fs::read_to_string(&config).unwrap();
            let dio_config: Config = toml::from_str(&contents).unwrap();
            let start_time = std::time::Instant::now();
            let params = DCRTPolyParams::new(
                dio_config.ring_dimension,
                dio_config.crt_depth,
                dio_config.crt_bits,
                dio_config.base_bits,
            );
            let log_base_q = params.modulus_digits();
            let switched_modulus = Arc::new(dio_config.switched_modulus);
            let dim = params.ring_dimension() as usize;
            assert_eq!(add_num + mul_num, log_base_q * 2);
            let public_circuit = BenchCircuit::new_add_mul(add_num, mul_num).as_poly_circuit();
            let public_circuit_num_input = public_circuit.num_input();
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

            let packed_input_size = obf_params.input_size.div_ceil(dim) + 1;
            assert_eq!(public_circuit_num_input, (2 * log_base_q) + packed_input_size - 1);
            let sampler_uniform = DCRTPolyUniformSampler::new();
            let mut rng = rand::rng();
            let hardcoded_key = sampler_uniform.sample_poly(&params, &DistType::BitDist);
            let obfuscation = obfuscate::<
                DCRTPolyMatrix,
                DCRTPolyUniformSampler,
                DCRTPolyHashSampler<Keccak256>,
                DCRTPolyTrapdoorSampler,
                _,
            >(obf_params.clone(), hardcoded_key.clone(), &mut rng);
            let obfuscation_time = start_time.elapsed();
            info!("Time to obfuscate: {:?}", obfuscation_time);
            let input = dio_config.input;
            assert_eq!(input.len(), dio_config.input_size);
            let output = obfuscation
                .eval::<DCRTPolyHashSampler<Keccak256>, DCRTPolyTrapdoorSampler>(
                    obf_params, &input,
                );
            let total_time = start_time.elapsed();

            info!("Time for evaluation: {:?}", total_time - obfuscation_time);
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

                let mut bool_eval = vec![];
                for e in eval {
                    let decompose_poly = e.to_bool_vec();
                    bool_eval.extend(decompose_poly);
                }

                assert_eq!(output, bool_eval);
            }
        }
    }
}
