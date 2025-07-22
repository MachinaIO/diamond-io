//! Public Lookup

use crate::{
    io::obf::store_and_drop_matrix,
    poly::{
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{collections::HashMap, path::Path};
use tokio::task::JoinHandle;
use tracing::info;

// const TAG_R_K: &[u8] = b"TAG_R_K";
// const TAG_A_PLT: &[u8] = b"A_PLT:";

/// Public Lookup Table
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PublicLut<P: Poly> {
    // // public matrix R_k that concated for all k
    // pub r_k_s: M,
    // public matrix A_z
    // a_z: Option<M>,
    // pub d: usize,
    /// m := (n+1)[logq]
    // m: usize,
    /// mapping f: k => (x_k, y_k)
    pub f: HashMap<usize, (P, P)>,
    //  Common public matrix A_LT: (n+1)xm
    // pub a_lt: M,
}

impl<P: Poly> PublicLut<P> {
    pub fn new(
        // params: &<M::P as Poly>::Params,
        // d: usize,
        f: HashMap<usize, (P, P)>,
        // r_k_hashkey: [u8; 32],
    ) -> Self {
        // let m = (1 + d) * params.modulus_digits();
        // let hash_sampler = SH::new();
        // let t = f.len();
        // todo: R_k could be sampled from uniform if we decided to dump it in disk.
        // let r_k_s = hash_sampler.sample_hash(
        //     params,
        //     r_k_hashkey,
        //     TAG_R_K,
        //     d + 1,
        //     m * t,
        //     DistType::FinRingDist,
        // );

        // let a_lt = hash_sampler.sample_hash(
        //     params,
        //     r_k_hashkey,
        //     TAG_A_PLT,
        //     d + 1,
        //     m,
        //     DistType::FinRingDist,
        // );
        // info!("A_LT ({}, {})", a_lt.row_size(), a_lt.col_size());
        // Self { r_k_s, a_z: None, f, d, m, a_lt }
        Self { f }
    }

    // /// Insert A_z public matrix
    // pub fn insert_a_z(&mut self, a_z: &M) {
    //     self.a_z = Some(a_z.clone())
    // }

    /// Find the row k with the maximum coefficient in the second M::P (y_k) of f HashMap
    /// Returns (k, max_coefficient)
    pub fn max_output_row(&self) -> (usize, <P as Poly>::Elem) {
        assert!(!self.f.is_empty(), "f must contain at least one element");
        self.f
            .iter()
            .filter_map(|(&k, (_, y_k))| y_k.coeffs().iter().max().cloned().map(|coeff| (k, coeff)))
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
        let matrices = (0..t)
            .into_par_iter()
            .map(|k| {
                info!("Processing k: {k}");
                let (x_k, y_k) = self.f.get(&k).expect("missing f(k)");
                let r_k = uniform_sampler.sample_uniform(params, d + 1, m, DistType::FinRingDist);
                info!("Sampled r_k ({}, {})", r_k.row_size(), r_k.col_size());
                let r_k_decomposed = r_k.decompose();
                let target_k = (r_k.clone() * x_k) + a_lt -
                    &(M::gadget_matrix(params, d + 1) * y_k) -
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
