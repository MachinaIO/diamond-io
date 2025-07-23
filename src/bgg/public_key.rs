use super::circuit::Evaluable;
use crate::{
    bgg::lut::public_lut::PublicLut,
    poly::{Poly, PolyMatrix},
    utils::debug_mem,
};
use rayon::prelude::*;
use std::{
    ops::{Add, Mul, Sub},
    path::PathBuf,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BggPublicKey<M: PolyMatrix> {
    pub matrix: M,
    pub reveal_plaintext: bool,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn new(matrix: M, reveal_plaintext: bool) -> Self {
        Self { matrix, reveal_plaintext }
    }

    pub fn concat_matrix(&self, others: &[Self]) -> M {
        self.matrix.concat_columns(&others.par_iter().map(|x| &x.matrix).collect::<Vec<_>>()[..])
    }

    // /// Writes the public key with id to files under the given directory.
    // pub async fn write_to_files<P: AsRef<std::path::Path> + Send + Sync>(
    //     &self,
    //     dir_path: P,
    //     id: &str,
    // ) {
    //     self.matrix.write_to_files(dir_path, id).await;
    // }

    /// Reads a public of given rows and cols with id from files under the given directory.
    pub fn read_from_files<P: AsRef<std::path::Path> + Send + Sync>(
        params: &<M::P as Poly>::Params,
        nrow: usize,
        ncol: usize,
        dir_path: P,
        id: &str,
        reveal_plaintext: bool,
    ) -> Self {
        let matrix = M::read_from_files(params, nrow, ncol, dir_path, id);
        Self { matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Add for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix: self.matrix + &other.matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Sub for BggPublicKey<M> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self - &other
    }
}

impl<M: PolyMatrix> Sub<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn sub(self, other: &Self) -> Self {
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix: self.matrix - &other.matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Mul for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        debug_mem(format!("BGGPublicKey::mul {:?}, {:?}", self.matrix.size(), other.matrix.size()));
        let decomposed = other.matrix.decompose();
        debug_mem("BGGPublicKey::mul decomposed");
        let matrix = self.matrix * decomposed;
        debug_mem("BGGPublicKey::mul matrix multiplied");
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        debug_mem("BGGPublicKey::mul reveal_plaintext");
        Self { matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Evaluable for BggPublicKey<M> {
    type Params = <M::P as Poly>::Params;
    type Matrix = M;

    fn rotate(self, params: &Self::Params, shift: usize) -> Self {
        debug_mem(format!("BGGPublicKey::rotate {:?}, {:?}", self.matrix.size(), shift));
        let rotate_poly = <M::P>::const_rotate_poly(params, shift);
        debug_mem("BGGPublicKey::rotate rotate_poly");
        let matrix = self.matrix.clone() * rotate_poly;
        debug_mem("BGGPublicKey::rotate matrix multiplied");
        Self { matrix, reveal_plaintext: self.reveal_plaintext }
    }

    fn from_digits(params: &Self::Params, one: &Self, digits: &[u32]) -> Self {
        debug_mem(format!("BGGPublicKey::from_digits {:?}, {:?}", one.matrix.size(), digits.len()));
        let const_poly =
            <M::P as Evaluable>::from_digits(params, &<M::P>::const_one(params), digits);
        debug_mem("BGGPublicKey::from_digits const_poly");
        let matrix = one.matrix.clone() * const_poly;
        debug_mem("BGGPublicKey::from_digits matrix multiplied");
        Self { matrix, reveal_plaintext: one.reveal_plaintext }
    }

    fn public_lookup(
        self,
        _: &Self::Params,
        plt: &mut PublicLut<Self::Matrix>,
        _: Option<(Self::Matrix, PathBuf, usize, usize)>,
    ) -> Self {
        let a_z = self.matrix;
        plt.insert_a_z(&a_z);
        Self { matrix: plt.a_lt.clone(), reveal_plaintext: self.reveal_plaintext }
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        bgg::{circuit::PolyCircuit, sampler::BGGPublicKeySampler},
        poly::dcrt::{params::DCRTPolyParams, DCRTPolyHashSampler},
        test_utils::setup_constant_plt,
    };
    use keccak_asm::Keccak256;

    #[test]
    fn test_pubkey_plt() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(1);
        let plt = setup_constant_plt(8, &params, d);
        let a_lt = plt.a_lt.clone();
        let plt_id = circuit.register_public_lookup(plt);
        let plt_gate = circuit.public_lookup_gate(inputs[0], plt_id);
        circuit.output(vec![plt_gate]);

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1], None);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, a_lt);
    }
}
