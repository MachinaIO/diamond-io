pub mod eval;
pub mod obf;
pub mod utils;

use crate::{bgg::BggEncoding, poly::matrix::*};

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix> {
    pub hash_key: [u8; 32],
    // pub fhe_enc: PolyMatrix<T, P, M>,
    pub b_fhe: M,
    pub encode_input: Vec<BggEncoding<M>>,
    pub encode_fhe_key: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages: Vec<(M, M)>,
    pub n_preimages: Vec<(M, M)>,
    pub k_preimages: Vec<(M, M)>,
    pub final_preimage: M,
}
