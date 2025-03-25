pub mod eval;
pub mod obf;
pub mod params;
pub mod utils;
use crate::{bgg::BggEncoding, poly::PolyMatrix};

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix> {
    pub hash_key: [u8; 32],
    pub enc_hardcoded_key: M,
    pub encodings_init: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages_ids: Vec<Vec<String>>,
    pub n_preimages_ids: Vec<Vec<String>>,
    pub k_preimages_ids: Vec<Vec<String>>,
    pub final_preimage_id: String,
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
