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
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PublicLut<M: PolyMatrix> {
    // public matrix R_k that concated for all k
    pub r_k_s: M,
    // public matrix A_z
    a_z: Option<M>,
    pub d: usize,
    /// m := (n+1)[logq]
    m: usize,
    /// mapping f: k => (x_k, y_k)
    pub f: HashMap<usize, (M::P, M::P)>,
    ///  Common public matrix A_LT: (n+1)xm
    pub a_lt: M,
}

impl<M: PolyMatrix> PublicLut<M> {
    pub fn new<SU: PolyUniformSampler<M = M>, SH: PolyHashSampler<[u8; 32], M = M>>(
        params: &<M::P as Poly>::Params,
        d: usize,
        f: HashMap<usize, (M::P, M::P)>,
        r_k_hashkey: [u8; 32],
    ) -> Self {
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

        let a_lt = hash_sampler.sample_hash(
            params,
            r_k_hashkey,
            TAG_A_PLT,
            d + 1,
            m,
            DistType::FinRingDist,
        );
        info!("A_LT ({}, {})", a_lt.row_size(), a_lt.col_size());
        Self { r_k_s, a_z: None, f, d, m, a_lt }
    }

    /// Insert A_z public matrix
    pub fn insert_a_z(&mut self, a_z: &M) {
        self.a_z = Some(a_z.clone())
    }

    /// Compute target, sample preimage and store it as file.
    pub fn preimage<ST>(
        &self,
        params: &<M::P as Poly>::Params,
        b_l: &M,
        b_l_plus_one: &M,
        trap_sampler: &ST,
        b_l_trapdoor: &ST::Trapdoor,
        b_l_plus_one_trapdoor: &ST::Trapdoor,
        input_size: usize,
        dir_path: &Path,
    ) -> Vec<JoinHandle<()>>
    where
        ST: PolyTrapdoorSampler<M = M> + Send + Sync,
        M: PolyMatrix + Send + 'static,
    {
        let t = self.f.len();
        let mut handles = Vec::new();

        let target_tuple: Vec<(usize, M)> = (0..t)
            .into_par_iter()
            .map(|k| {
                let (x_k, y_k) = self.f.get(&k).expect("missing f(k)");
                let r_k = self.r_k_s.slice_columns(k * self.m, (k + 1) * self.m);
                // rhs  = A_LT - A_z·G⁻¹(R_k) + x_k·R_k - y_k·G
                let rhs = self.a_lt.clone() + (r_k.clone() * x_k) -
                    &(M::gadget_matrix(params, self.d + 1) * y_k) -
                    self.a_z.clone().expect("A_z should be exist") * r_k.decompose();
                (k, rhs)
            })
            .collect();

        // first sample L_common
        let id = M::identity(params, b_l_plus_one.row_size(), None);
        let zeros = M::zero(params, (input_size - 1) * id.row_size(), id.col_size());
        let tensor_lhs = id.concat_rows(&[&zeros]);
        // info!("id ({}, {})", id.row_size(), id.col_size());
        // info!("tensor_lhs ({}, {})", tensor_lhs.row_size(), tensor_lhs.col_size());
        // info!("b_l_plus_one ({}, {})", b_l_plus_one.row_size(), b_l_plus_one.col_size());
        // info!("b_l ({}, {})", b_l.row_size(), b_l.col_size());
        let l_common_target = tensor_lhs * b_l_plus_one;
        let l_common = trap_sampler.preimage(params, b_l_trapdoor, b_l, &l_common_target);
        handles.push(store_and_drop_matrix(l_common, dir_path, &format!("L_common")));

        for (k, target_k) in target_tuple {
            // let zeros = M::zero(params, (input_size - 1) * rhs_k.row_size(), rhs_k.col_size());
            // let target = rhs_k.concat_rows(&[&zeros]);
            info!("target_k ({}, {})", target_k.row_size(), target_k.col_size());
            let l_k = trap_sampler.preimage(params, b_l_plus_one_trapdoor, b_l_plus_one, &target_k);
            handles.push(store_and_drop_matrix(l_k, dir_path, &format!("L_{k}")));
        }

        handles
    }
}
