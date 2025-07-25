use super::circuit::Evaluable;
use crate::{
    bgg::{circuit::PltEvaluator, lut::public_lut::PublicLut},
    poly::{
        sampler::{PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix,
    },
    utils::debug_mem,
};
use rayon::prelude::*;
use std::{
    marker::PhantomData,
    ops::{Add, Mul, Sub},
    path::PathBuf,
    sync::Arc,
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
    type P = M::P;

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
}

#[derive(Debug)]
pub struct BggPubKeyPltEvaluator<M, SH, SU, ST>
where
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Send + Sync,
{
    pub hash_key: [u8; 32],
    pub trap_sampler: ST,
    pub pub_matrix: Arc<M>,
    pub trapdoor: Arc<ST::Trapdoor>,
    pub dir_path: PathBuf,
    _sh: PhantomData<SH>,
    _su: PhantomData<SU>,
    _st: PhantomData<ST>,
}

impl<M, SH, SU, ST> PltEvaluator<BggPublicKey<M>> for BggPubKeyPltEvaluator<M, SH, SU, ST>
where
    M: PolyMatrix + Send + 'static,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Send + Sync,
{
    fn public_lookup(
        &self,
        params: &<BggPublicKey<M> as Evaluable>::Params,
        plt: &PublicLut<<BggPublicKey<M> as Evaluable>::P>,
        input: BggPublicKey<M>,
        id: usize,
    ) -> BggPublicKey<M> {
        let d = input.matrix.row_size() - 1;
        let a_lt = plt.derive_a_lt::<M, SH>(params, d, self.hash_key, id);
        plt.preimage::<M, SU, ST>(
            params,
            &self.trap_sampler,
            &self.pub_matrix,
            &self.trapdoor,
            &input.matrix,
            &a_lt,
            id,
            &self.dir_path,
        );
        BggPublicKey { matrix: a_lt, reveal_plaintext: true }
    }
}

impl<M, SH, SU, ST> BggPubKeyPltEvaluator<M, SH, SU, ST>
where
    M: PolyMatrix,
    SH: PolyHashSampler<[u8; 32], M = M> + Send + Sync,
    SU: PolyUniformSampler<M = M> + Send + Sync,
    ST: PolyTrapdoorSampler<M = M> + Send + Sync,
{
    pub fn new(
        hash_key: [u8; 32],
        trap_sampler: ST,
        pub_matrix: Arc<M>,
        trapdoor: Arc<ST::Trapdoor>,
        dir_path: PathBuf,
    ) -> Self {
        Self {
            hash_key,
            trap_sampler,
            pub_matrix,
            trapdoor,
            dir_path,
            _sh: PhantomData,
            _su: PhantomData,
            _st: PhantomData,
        }
    }
}
