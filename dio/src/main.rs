use clap::Parser;
use diamond_io::{
    bgg::circuit::serde::SerializablePolyCircuit,
    io::{obf::obfuscate, params::ObfuscationParams, serde::SerializableObfuscationParams},
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
use rand::Rng;
use serde::Deserialize;
use std::{
    fs::File,
    io::{BufReader, Read},
    path::PathBuf,
    sync::Arc,
};
use tracing::info;

/// Simple program to obfuscate and evaluate
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    config: PathBuf,
}

#[derive(Clone, Debug, Deserialize)]
struct Config {
    dcrt_poly_params: DCRTPolyParams,
    obfuscation_params: SerializableObfuscationParams,
}

fn main() {
    init_tracing();
    let args = Args::parse();
    let file = File::open(args.config).expect("invalid path");
    let mut buf_reader = BufReader::new(file);
    let mut contents = String::new();
    buf_reader.read_to_string(&mut contents).unwrap();
    let dio_config: Config = toml::from_str(&contents).unwrap();
    let start_time = std::time::Instant::now();
    let params = DCRTPolyParams::new(
        dio_config.dcrt_poly_params.ring_dimension(),
        dio_config.dcrt_poly_params.crt_depth(),
        dio_config.dcrt_poly_params.crt_bits(),
        dio_config.dcrt_poly_params.base_bits(),
    );

    let switched_modulus = Arc::new(dio_config.obfuscation_params.switched_modulus);
    let file = File::open(dio_config.obfuscation_params.public_circuit_path).unwrap();
    let mut buf_reader = BufReader::new(file);
    let mut contents = String::new();
    buf_reader.read_to_string(&mut contents).unwrap();
    let serialized_public_circuit: SerializablePolyCircuit =
        serde_json::from_str(&contents).expect("JSON was not well-formatted");
    let public_circuit = serialized_public_circuit.to_circuit();

    let obf_params = ObfuscationParams {
        params: params.clone(),
        switched_modulus,
        input_size: dio_config.obfuscation_params.input_size,
        level_width: dio_config.obfuscation_params.level_width,
        public_circuit: public_circuit.clone(),
        d: dio_config.obfuscation_params.d,
        encoding_sigma: dio_config.obfuscation_params.encoding_sigma,
        hardcoded_key_sigma: dio_config.obfuscation_params.hardcoded_key_sigma,
        p_sigma: dio_config.obfuscation_params.p_sigma,
        trapdoor_sigma: dio_config.obfuscation_params.trapdoor_sigma,
    };

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

    let bool_in = rng.random::<bool>();
    let input = [bool_in, false, false, false];
    let output = obfuscation
        .eval::<DCRTPolyHashSampler<Keccak256>, DCRTPolyTrapdoorSampler>(obf_params, &input);
    let total_time = start_time.elapsed();
    info!("Time for evaluation: {:?}", total_time - obfuscation_time);
    info!("Total time: {:?}", total_time);
    let n = output.len() / 2;
    let output_1st_gate = output[..n].to_vec();
    let output_2nd_gate = output[n..].to_vec();
    let input_poly =
        DCRTPoly::from_const(&params, &FinRingElem::constant(&params.modulus(), bool_in as u64));
    assert_eq!(output_1st_gate, (hardcoded_key.clone() * input_poly.clone()).to_bool_vec());
    assert_eq!(output_2nd_gate, (hardcoded_key * input_poly).to_bool_vec());
}
