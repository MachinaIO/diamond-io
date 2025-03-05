// pub mod ciphertext;
// pub mod circuit;
// pub mod operations;
// pub mod parameters;
pub mod circuits;
pub mod eval;
pub mod sampler;
use crate::poly::matrix::*;
use itertools::Itertools;
// use thiserror::Error;

// #[derive(Error, Debug)]
// pub enum BggError {
//     #[error("Unknown plaintext for the left-hand input of multiplication: left-input index {0}, output index {1}")]
//     UnknownPlaintextForMul(usize, usize),
// }

#[derive(Debug, Clone)]
pub struct BggPublicKey<M: PolyMatrix> {
    pub matrix: M,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn new(matrix: M) -> Self {
        Self { matrix }
    }

    pub fn concat_matrix(&self, others: &[Self]) -> M {
        self.matrix.concat_columns(&others.iter().map(|x| x.matrix.clone()).collect_vec()[..])
    }
}

#[derive(Debug, Clone)]
pub struct BggEncoding<M: PolyMatrix> {
    pub vector: M,
    pub pubkey: BggPublicKey<M>,
    pub plaintext: Option<<M as PolyMatrix>::P>,
}

impl<M: PolyMatrix> BggEncoding<M> {
    pub fn new(
        vector: M,
        pubkey: BggPublicKey<M>,
        plaintext: Option<<M as PolyMatrix>::P>,
    ) -> Self {
        Self { vector, pubkey, plaintext }
    }

    pub fn concat_vector(&self, others: &[Self]) -> M {
        self.vector.concat_columns(&others.iter().map(|x| x.vector.clone()).collect_vec()[..])
    }
}
