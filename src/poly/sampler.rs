use std::path::PathBuf;

use super::{Poly, PolyMatrix};

#[derive(Debug)]
/// Enum representing different types of distributions for random sampling.
pub enum DistType {
    /// Distribution over a finite ring, typically samples elements from a ring in a uniform or
    /// near-uniform manner
    FinRingDist,
    /// discrete Gaussian distribution described in [[GPV08](https://eprint.iacr.org/2007/432),[BDJ+24](https://eprint.iacr.org/2024/1742)], where
    /// noise is drawn from a discrete Gaussian over a lattice Λ with parameter σ > 0.
    /// Each sample is drawn proportionally to exp(-π‖x‖² / σ²), restricted to x ∈ Λ.
    ///
    /// * `sigma` - The Gaussian parameter (standard deviation).
    GaussDist { sigma: f64 },
    /// Distribution that produces random bits (0 or 1).
    BitDist,
}

/// Trait for sampling a polynomial based on a hash function.
pub trait PolyHashSampler<K: AsRef<[u8]>> {
    type M: PolyMatrix;
    /// Samples a matrix of ring elements from a pseudorandom source defined by a hash function `H`
    /// Compute H(key || tag || i)
    ///
    /// and a distribution type specified by `dist`.
    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        tag: B,
        nrow: usize,
        ncol: usize,
        dist: DistType,
    ) -> Self::M;

    // Set the key of the sampler.
    fn set_key(&mut self, key: K);
}

pub trait PolyUniformSampler {
    type M: PolyMatrix;

    fn sample_uniform(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        nrow: usize,
        ncol: usize,
        dist: DistType,
    ) -> Self::M;
}

pub trait PolyTrapdoorSampler {
    type M: PolyMatrix;
    type Trapdoor;
    type MatrixPtr;

    fn trapdoor(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        size: usize,
    ) -> (Self::Trapdoor, Self::M);
    fn preimage_to_fs(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
        preimage_id: &str,
    ) -> Vec<PathBuf>;
    fn process_preimage_block_to_fs(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::MatrixPtr,
        target_block: &Self::M,
        size: usize,
        preimage_block_id: &str,
    ) -> PathBuf;
    fn preimage_from_fs(
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        preimages_paths: &[PathBuf],
    ) -> Self::M;
}
