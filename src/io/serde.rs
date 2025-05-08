use std::{path::PathBuf, str::FromStr};

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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableObfuscationParams {
    #[serde(serialize_with = "biguint_to_string", deserialize_with = "biguint_from_string")]
    pub switched_modulus: BigUint,
    pub input_size: usize,
    pub level_width: usize,
    pub public_circuit_path: PathBuf,
    pub d: usize,
    pub hardcoded_key_sigma: f64,
    pub p_sigma: f64,
    pub trapdoor_sigma: f64,
}
