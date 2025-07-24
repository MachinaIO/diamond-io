//! Public Lookup

use crate::{
    io::obf::store_and_drop_matrix,
    poly::{
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use rayon::prelude::*;
use std::{collections::HashMap, path::Path};
use tokio::task::JoinHandle;
use tracing::info;

/// Public Lookup Table
#[derive(Debug, Clone, Default)]
pub struct PublicLut<P: Poly> {
    pub f: HashMap<P, (usize, P)>,
}

impl<P: Poly> PublicLut<P> {
    pub fn new(f: HashMap<P, (usize, P)>) -> Self {
        Self { f }
    }

    /// Find the row k with the maximum coefficient in the second M::P (y_k) of f HashMap
    /// Returns (k, max_coefficient)
    pub fn max_output_row(&self) -> (usize, <P as Poly>::Elem) {
        assert!(!self.f.is_empty(), "f must contain at least one element");
        self.f
            .iter()
            .filter_map(|(_, (k, y_k))| y_k.coeffs().iter().max().cloned().map(|coeff| (*k, coeff)))
            .max_by(|a, b| a.1.cmp(&b.1))
            .expect("no coefficients found in any y_k")
    }

    pub fn derive_a_lt<M, SH>(
        &self,
        params: &<M::P as Poly>::Params,
        d: usize,
        hash_key: [u8; 32],
        id: usize,
    ) -> M
    where
        M: PolyMatrix<P = P>,
        SH: PolyHashSampler<[u8; 32], M = M>,
    {
        info!("Deriving A_LT for id: {id}");
        let m = (d + 1) * params.modulus_digits();
        let hash_sampler = SH::new();
        let tag = format!("A_LT_{id}");
        info!("Tag for A_LT: {tag}");
        hash_sampler.sample_hash(
            params,
            hash_key,
            tag.into_bytes(),
            d + 1,
            m,
            DistType::FinRingDist,
        )
    }

    /// Compute target, sample preimage and store it as file.
    pub fn preimage<M, SU, ST>(
        &self,
        params: &<M::P as Poly>::Params,
        trap_sampler: &ST,
        pub_matrix: &M,
        trapdoor: &ST::Trapdoor,
        a_z: &M,
        a_lt: &M,
        id: usize,
        dir_path: &Path,
        handles_out: &mut Vec<JoinHandle<()>>,
    ) where
        M: PolyMatrix<P = P> + Send + 'static,
        SU: PolyUniformSampler<M = M> + Send + Sync,
        ST: PolyTrapdoorSampler<M = M> + Send + Sync,
    {
        info!("Preimage for id: {id}");
        let t = self.f.len();
        let d = pub_matrix.row_size() - 1;
        let m = (d + 1) * params.modulus_digits();
        let uniform_sampler = SU::new();
        let matrices = self
            .f
            .iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(x_k, (k, y_k))| {
                info!("Processing k: {k}");
                let r_k = uniform_sampler.sample_uniform(params, d + 1, m, DistType::FinRingDist);
                info!("Sampled r_k ({}, {})", r_k.row_size(), r_k.col_size());
                let r_k_decomposed = r_k.decompose();
                let target_k = (r_k.clone() * x_k) + a_lt -
                    &(pub_matrix.clone() * y_k) -
                    a_z.clone() * r_k_decomposed;
                info!("target_k ({}, {})", target_k.row_size(), target_k.col_size());
                (r_k, trap_sampler.preimage(params, trapdoor, pub_matrix, &target_k))
            })
            .collect::<Vec<_>>();
        info!("Preimage matrices computed for id: {id}");
        // [TODO] Use a channel within the above iterator to bound the memory usage.
        for (k, (r_k, l_k)) in matrices.into_iter().enumerate() {
            handles_out.push(store_and_drop_matrix(r_k, dir_path, &format!("R_{id}_{k}")));
            handles_out.push(store_and_drop_matrix(l_k, dir_path, &format!("L_{id}_{k}")));
        }
    }
}
