use std::path::PathBuf;

use num_bigint::BigUint;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableObfuscationParams {
    pub switched_modulus: BigUint,
    pub input_size: usize,
    pub level_width: usize,
    pub public_circuit_path: PathBuf,
    pub d: usize,
    pub encoding_sigma: f64,
    pub hardcoded_key_sigma: f64,
    pub p_sigma: f64,
    pub trapdoor_sigma: f64,
}
