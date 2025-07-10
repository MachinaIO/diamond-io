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

const TAG_R_K: &[u8] = b"TAG_R_K";
const TAG_A_PLT: &[u8] = b"A_PLT:";

/// Public Lookup Table
/// Considering adjusting on diamond-io case.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PublicLut<M: PolyMatrix> {
    // public matrix different for all k: (n+1)xm
    pub r_k_s: M,
    pub target_hashmap: HashMap<usize, M>,
    d: usize,
    m: usize,
    // todo: yes i know potentially we should sample r_k via hash sampler and slice thru k but for
    // k => (x_k, y_k)
    pub f: HashMap<usize, (M::P, M::P)>,
    pub a_lt: M,
}

impl<M: PolyMatrix> PublicLut<M> {
    pub fn new<SU: PolyUniformSampler<M = M>, SH: PolyHashSampler<[u8; 32], M = M>>(
        params: &<M::P as Poly>::Params,
        d: usize,
        f: HashMap<usize, (M::P, M::P)>,
        r_k_hashkey: [u8; 32],
    ) -> Self {
        // m := (n+1)[logq]
        let m = (1 + d) * params.modulus_digits();
        let hash_sampler = SH::new();
        let t = f.len();
        let r_k_s = hash_sampler.sample_hash(
            params,
            r_k_hashkey,
            TAG_R_K,
            d + 1,
            m * t,
            DistType::FinRingDist,
        );

        // new public BGG+matrix common for all rows: (n+1)xm
        let a_lt = hash_sampler.sample_hash(
            params,
            r_k_hashkey,
            TAG_A_PLT,
            d + 1,
            m,
            DistType::FinRingDist,
        );
        info!("A_LT ({}, {})", a_lt.row_size(), a_lt.col_size());
        Self { r_k_s, target_hashmap: Default::default(), f, d, m, a_lt }
    }

    /// interface will be called in eval for compute rhs
    pub fn compute_target(&mut self, params: &<M::P as Poly>::Params, a_z: &M) {
        let t = self.f.len();

        let target_tuple: Vec<(usize, M)> = (0..t)
            .into_par_iter()
            .map(|k| {
                let (x_k, y_k) = self.f.get(&k).expect("missing f(k)");
                let r_k = self.r_k_s.slice_columns(k * self.m, (k + 1) * self.m);
                // rhs  = A_LT - A_z·G⁻¹(R_k) + x_k·R_k - y_k·G
                let rhs = self.a_lt.clone() + (r_k.clone() * x_k) -
                    &(M::gadget_matrix(params, self.d + 1) * y_k) -
                    a_z.clone() * r_k.decompose();
                (k, rhs)
            })
            .collect();
        self.target_hashmap = target_tuple.into_iter().collect();
    }

    /// interface will be called in diamond io for storing preimage
    pub fn preimage<ST>(
        &self,
        params: &<M::P as Poly>::Params,
        b_l: &M,
        trap_sampler: &ST,
        trapdoor: &ST::Trapdoor,
        input_size: usize,
        dir_path: &Path,
    ) -> Vec<JoinHandle<()>>
    where
        ST: PolyTrapdoorSampler<M = M> + Send + Sync,
        M: PolyMatrix + Send + 'static,
    {
        assert_eq!(self.f.len(), self.target_hashmap.len());

        let mut handles = Vec::new();
        for (k, rhs_k) in &self.target_hashmap {
            let zeros = M::zero(params, (input_size - 1) * rhs_k.row_size(), rhs_k.col_size());
            let target = rhs_k.concat_rows(&[&zeros]);
            let l_k = trap_sampler.preimage(params, trapdoor, b_l, &target);
            handles.push(store_and_drop_matrix(l_k, dir_path, &format!("L_{k}")));
        }
        handles
    }
}
