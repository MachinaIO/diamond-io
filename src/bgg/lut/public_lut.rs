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

// #[cfg(test)]
// mod tests {
//     use keccak_asm::Keccak256;

//     use super::*;
//     use crate::{
//         bgg::lut::utils::p_vector_for_inputs,
//         poly::{
//             dcrt::{
//                 DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams,
//                 DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler,
//             },
//             sampler::PolyTrapdoorSampler,
//         },
//         utils::init_tracing,
//     };

//     #[test]
//     fn test_public_lut() {
//         init_tracing();

//         /* Setup */
//         let params = DCRTPolyParams::default();
//         let uni = DCRTPolyUniformSampler::new();
//         let d = 1;
//         // let input_size = 2;

//         /* Lookup mapping k => (x_k, y_k) */
//         let mut f = HashMap::new();
//         f.insert(0, (DCRTPoly::const_int(&params, 0), DCRTPoly::const_int(&params, 7)));
//         f.insert(1, (DCRTPoly::const_int(&params, 1), DCRTPoly::const_int(&params, 5)));
//         f.insert(2, (DCRTPoly::const_int(&params, 2), DCRTPoly::const_int(&params, 6)));
//         f.insert(3, (DCRTPoly::const_int(&params, 3), DCRTPoly::const_int(&params, 1)));
//         f.insert(4, (DCRTPoly::const_int(&params, 4), DCRTPoly::const_int(&params, 0)));
//         f.insert(5, (DCRTPoly::const_int(&params, 5), DCRTPoly::const_int(&params, 3)));
//         f.insert(6, (DCRTPoly::const_int(&params, 6), DCRTPoly::const_int(&params, 4)));
//         f.insert(7, (DCRTPoly::const_int(&params, 7), DCRTPoly::const_int(&params, 2)));

//         /* Obfuscation Step */
//         // let trap_sampler = DCRTPolyTrapdoorSampler::new(&params, 4.578);
//         // let (trapdoor, b_l) = trap_sampler.trapdoor(&params, (1 + input_size) * (d + 1));
//         let t = f.len();
//         /* BGG+ encoding setup */
//         let secrets = uni.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
//         // in reality there should be input insertion step that updates the secret
//         let s_x_l = {
//             let minus_one_poly = DCRTPoly::const_minus_one(&params);
//             let mut secrets = secrets.to_vec();
//             secrets.push(minus_one_poly);
//             DCRTPolyMatrix::from_poly_vec_row(&params, secrets)
//         };
//         info!("s_x_l ({},{})", s_x_l.row_size(), s_x_l.col_size());

//         let lut = PublicLut::<DCRTPolyMatrix>::new::<
//             DCRTPolyUniformSampler,
//             DCRTPolyHashSampler<Keccak256>,
//         >(&params, d, f);

//         // /* Evaluation Step */
//         // let inputs = vec![1, 0];
//         // let k = 6;
//         // assert_eq!(inputs.len(), input_size);

//         // let p_x_l = p_vector_for_inputs(&b_l, inputs, &params, 0.0, &lut.s_x_l);
//         // let c_y_k = lut.evaluate(&params, t, p_x_l, k);
//         // let (_, y_k) = f.get(&k).expect("missing f(k)");
//         // let expected_c_y_k =
//         //     &lut.s_x_l * (&lut.a_lt - &(DCRTPolyMatrix::gadget_matrix(&params, d + 1) * y_k));
//         // debug_assert_eq!(expected_c_y_k, c_y_k);
//     }
// }
