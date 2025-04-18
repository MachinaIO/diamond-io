#[cfg(feature = "bgm")]
pub mod bgm;

use crate::{bgg::BggEncoding, poly::PolyMatrix};

pub mod eval;
pub mod obf;
pub mod params;
pub mod serde;
pub mod utils;

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix> {
    pub hash_key: [u8; 32],
    pub b: M,
    pub encodings_init: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages: Vec<Vec<M>>,
    pub n_preimages: Vec<Vec<M>>,
    pub k_preimages: Vec<Vec<M>>,
    pub final_preimage: M,
    #[cfg(feature = "debug")]
    pub s_init: M,
    #[cfg(feature = "debug")]
    pub minus_t_bar: <M as PolyMatrix>::P,
    #[cfg(feature = "debug")]
    pub bs: Vec<Vec<M>>,
    #[cfg(feature = "debug")]
    pub hardcoded_key: <M as PolyMatrix>::P,
}
