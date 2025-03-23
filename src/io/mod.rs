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
    pub m_preimages: Vec<(M, M)>,
    pub n_preimages: Vec<(M, M)>,
    pub k_preimages: Vec<(M, M)>,
    pub final_preimage: M,
    #[cfg(feature = "test")]
    pub s_init: M,
    #[cfg(feature = "test")]
    pub t_bar: <M as PolyMatrix>::P,
    #[cfg(feature = "test")]
    pub bs: Vec<(M, M, M)>,
    #[cfg(feature = "test")]
    pub hardcoded_key: <M as PolyMatrix>::P,
    #[cfg(feature = "test")]
    pub final_preimage_target: M,
}
