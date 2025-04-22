use std::str::FromStr;

use num_bigint::BigUint;
use serde::{de, Deserialize, Deserializer, Serialize, Serializer};

fn biguint_to_string<S>(value: &BigUint, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&value.to_str_radix(10))
}

fn biguint_from_string<'de, D>(deserializer: D) -> Result<BigUint, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    BigUint::from_str(&s).map_err(de::Error::custom)
}

fn default_trapdoor_sigma() -> Option<f64> {
    Some(4.578)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    #[serde(serialize_with = "biguint_to_string", deserialize_with = "biguint_from_string")]
    pub switched_modulus: BigUint,
    pub input_size: usize,
    pub level_width: usize,
    pub d: usize,
    pub encoding_sigma: f64,
    pub hardcoded_key_sigma: f64,
    pub p_sigma: f64,
    #[serde(default = "default_trapdoor_sigma")]
    pub trapdoor_sigma: Option<f64>,
    /// polynomial ring dimension
    pub ring_dimension: u32,
    /// size of the tower
    pub crt_depth: usize,
    /// number of bits of each tower's modulus
    pub crt_bits: usize,
    /// bit size of the base for the gadget vector and decomposition
    pub base_bits: u32,
    pub input: Vec<bool>,
}
