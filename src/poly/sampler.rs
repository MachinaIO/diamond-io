use super::{Poly, PolyMatrix};

#[derive(Debug, Clone, Copy)]
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

    fn new() -> Self;

    /// Samples a matrix of ring elements from a pseudorandom source defined by a hash function `H`
    /// Compute H(key || tag || i)
    ///
    /// and a distribution type specified by `dist`.
    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        key: [u8; 32],
        tag: B,
        nrow: usize,
        ncol: usize,
        dist: DistType,
    ) -> Self::M;
}

pub trait PolyUniformSampler {
    type M: PolyMatrix;

    fn new() -> Self;

    fn sample_poly(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        dist: &DistType,
    ) -> <Self::M as PolyMatrix>::P;

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
    type Trapdoor: Send + Sync;

    fn new(params: &<<Self::M as PolyMatrix>::P as Poly>::Params, sigma: f64) -> Self;

    fn trapdoor(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        size: usize,
    ) -> (Self::Trapdoor, Self::M);

    fn preimage(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
    ) -> Self::M;
}
