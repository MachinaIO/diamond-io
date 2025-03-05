use super::{Poly, PolyMatrix};

#[derive(Debug)]
/// Enum representing different types of distributions for random sampling.
pub enum DistType {
    /// Distribution over a finite ring, typically samples elements from a ring in a uniform or near-uniform manner
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

pub trait PolySampler {
    type M: PolyMatrix + Send + Sync;

    fn get_params(&self) -> <<Self::M as PolyMatrix>::P as Poly>::Params;
}

/// Trait for sampling a polynomial based on a hash function.
pub trait PolyHashSampler<K: AsRef<[u8]>>: PolySampler {
    /// Samples a matrix of ring elements from a pseudorandom source defined by a hash function `H` Compute H(key || tag || i)
    ///
    /// and a distribution type specified by `dist`.
    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        tag: B,
        rows: usize,
        columns: usize,
        dist: DistType,
    ) -> Self::M;

    // Set the key of the sampler.
    fn set_key(&mut self, key: K);
}

pub trait PolyUniformSampler: PolySampler {
    fn sample_uniform(&self, rows: usize, columns: usize, dist: DistType) -> Self::M;
}

pub trait PolyTrapdoorSampler: PolySampler {
    type Trapdoor;
    fn trapdoor(&self) -> (Self::Trapdoor, Self::M);
    fn preimage(
        &self,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
    ) -> Self::M;
}
