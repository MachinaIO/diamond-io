#[cfg(feature = "bgm")]
pub mod bgm;
pub mod eval;
pub mod obf;
pub mod params;
pub mod utils;
use std::path::PathBuf;

use crate::{bgg::BggEncoding, poly::PolyMatrix};

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix> {
    pub hash_key: [u8; 32],
    pub enc_hardcoded_key: M,
    pub encodings_init: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages_paths: Vec<Vec<Vec<PathBuf>>>,
    pub n_preimages_paths: Vec<Vec<Vec<PathBuf>>>,
    pub k_preimages_paths: Vec<Vec<Vec<PathBuf>>>,
    pub final_preimage_path: Vec<PathBuf>,
    #[cfg(feature = "test")]
    pub s_init: M,
    #[cfg(feature = "test")]
    pub t_bar: <M as PolyMatrix>::P,
    #[cfg(feature = "test")]
    pub bs: Vec<Vec<M>>,
    #[cfg(feature = "test")]
    pub hardcoded_key: <M as PolyMatrix>::P,
    #[cfg(feature = "test")]
    pub final_preimage_target: M,
}
