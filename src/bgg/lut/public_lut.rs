//! Public Lookup

use crate::{
    poly::{
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
    storage::store_and_drop_matrix,
    utils::log_mem,
};
use rayon::prelude::*;
use std::{collections::HashMap, path::Path};

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
        log_mem(format!("Deriving A_LT for id: {id}"));
        let m = (d + 1) * params.modulus_digits();
        let hash_sampler = SH::new();
        let tag = format!("A_LT_{id}");
        log_mem(format!("Tag for A_LT: {tag}"));
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
    ) where
        M: PolyMatrix<P = P> + Send + 'static,
        SU: PolyUniformSampler<M = M> + Send + Sync,
        ST: PolyTrapdoorSampler<M = M> + Send + Sync,
    {
        let d = pub_matrix.row_size() - 1;
        let m = (d + 1) * params.modulus_digits();
        let uniform_sampler = SU::new();
        let gadget = M::gadget_matrix(params, d + 1);
        let matrices = self
            .f
            .iter()
            .collect::<Vec<_>>()
            .into_par_iter()
            .map(|(x_k, (k, y_k))| {
                let r_k = uniform_sampler.sample_uniform(params, d + 1, m, DistType::FinRingDist);
                let target_k = (r_k.clone() * x_k) + a_lt -
                    &(gadget.clone() * y_k) -
                    (a_z.clone() * r_k.decompose());
                (k, r_k, trap_sampler.preimage(params, trapdoor, pub_matrix, &target_k))
            })
            .collect::<Vec<_>>();
        log_mem(format!("Preimage matrices computed for id: {id}"));
        for (k, r_k, l_k) in matrices.into_iter() {
            store_and_drop_matrix(r_k, dir_path, &format!("R_{id}_{k}"));
            store_and_drop_matrix(l_k, dir_path, &format!("L_{id}_{k}"));
        }
    }
}
