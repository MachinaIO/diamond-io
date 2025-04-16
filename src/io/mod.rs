#[cfg(feature = "bgm")]
pub mod bgm;

use crate::{bgg::BggEncoding, poly::PolyMatrix};

pub mod eval;
pub mod obf;
pub mod params;
pub mod utils;

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix, P: AsRef<std::path::Path> + Send + Sync = std::path::PathBuf>
{
    pub hash_key: [u8; 32],
    pub ct_b: M,
    pub encodings_init: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages: Vec<Vec<M>>,
    pub n_preimages: Vec<Vec<M>>,
    pub k_preimages: Vec<Vec<M>>,
    pub final_preimage: M,
    pub dir_path: P,
    #[cfg(feature = "test")]
    pub s_init: M,
    #[cfg(feature = "test")]
    pub minus_t_bar: <M as PolyMatrix>::P,
    #[cfg(feature = "test")]
    pub bs: Vec<Vec<M>>,
    #[cfg(feature = "test")]
    pub hardcoded_key: <M as PolyMatrix>::P,
    #[cfg(feature = "test")]
    pub final_preimage_target: M,
}
