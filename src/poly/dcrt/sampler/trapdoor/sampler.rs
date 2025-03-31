use rayon::iter::ParallelIterator;

use super::{
    trapdoor::DCRTTrapdoor,
    utils::{gen_dgg_int_vec, gen_int_karney, split_int64_vec_to_elems},
};
use crate::{
    parallel_iter,
    poly::{
        dcrt::{
            matrix::{i64_matrix::I64MatrixParams, I64Matrix},
            sampler::DCRTPolyUniformSampler,
            DCRTPoly, DCRTPolyMatrix, DCRTPolyParams,
        },
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        PolyMatrix, PolyParams,
    },
};
use std::path::PathBuf;

const SIGMA: f64 = 4.578;
const SPECTRAL_CONSTANT: f64 = 1.8;

pub struct DCRTPolyTrapdoorSampler {
    sigma: f64,
    c: f64,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(sigma: f64) -> Self {
        // base = 2
        let c = 3.0 * SIGMA;
        Self { sigma, c }
    }

    fn preimage_square(
        &self,
        params: &DCRTPolyParams,
        trapdoor: &DCRTTrapdoor,
        public_matrix: &DCRTPolyMatrix,
        target: &DCRTPolyMatrix,
        s: f64,
        dgg_large: (f64, f64),
        peikert: bool,
    ) -> DCRTPolyMatrix {
        // (d * (k+2)) times d
        let p_hat = trapdoor.sample_pert_square_mat(s, self.c, self.sigma, dgg_large, peikert);

        let perturbed_syndrome = target.clone() - public_matrix.clone() * &p_hat;
        let k = params.modulus_bits();
        let d = public_matrix.row_size();
        let depth = params.crt_depth();

        let mut z_hat_vecs = Vec::with_capacity(d);
        for i in 0..d {
            let mut z_hat_row_vec = Vec::with_capacity(d);
            for j in 0..d {
                let z_hat_bbi_blocks = parallel_iter!(0..depth)
                    .map(|tower_idx| {
                        gauss_samp_gq_arb_base(
                            params,
                            tower_idx,
                            &perturbed_syndrome.entry(i, j),
                            self.c,
                            self.sigma,
                        )
                    })
                    .collect::<Vec<_>>();
                let z_hat_bbi = z_hat_bbi_blocks[0]
                    .concat_rows(&z_hat_bbi_blocks[1..].iter().collect::<Vec<_>>());
                let z_hat = split_int64_vec_to_elems(&z_hat_bbi, params);
                z_hat_row_vec.push(z_hat);
            }
            z_hat_vecs.push(
                z_hat_row_vec[0].concat_columns(&z_hat_row_vec[1..].iter().collect::<Vec<_>>()),
            );
        }
        let z_hat_mat = z_hat_vecs[0].concat_rows(&z_hat_vecs[1..].iter().collect::<Vec<_>>());

        let r_z_hat = trapdoor.r.clone() * &z_hat_mat;
        let e_z_hat = trapdoor.e.clone() * &z_hat_mat;
        let p_hat_former = (p_hat.slice_rows(0, d) + r_z_hat)
            .concat_rows(&[&(p_hat.slice_rows(d, 2 * d) + e_z_hat)]);
        let p_hat_latter = p_hat.slice_rows(2 * d, d * (k + 2)) + z_hat_mat;
        p_hat_former.concat_rows(&[&p_hat_latter])
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = DCRTTrapdoor;

    fn trapdoor(
        &self,
        params: &<<Self::M as crate::poly::PolyMatrix>::P as crate::poly::Poly>::Params,
        size: usize,
    ) -> (Self::Trapdoor, Self::M) {
        let trapdoor = DCRTTrapdoor::new(params, size, self.sigma);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let a_bar = uniform_sampler.sample_uniform(params, size, size, DistType::FinRingDist);
        let g = DCRTPolyMatrix::gadget_matrix(params, size);
        let a0 = a_bar.concat_columns(&[&DCRTPolyMatrix::identity(params, size, None)]);
        let a1 = g - (a_bar * &trapdoor.r + &trapdoor.e);
        let a = a0.concat_columns(&[&a1]);
        (trapdoor, a)
    }

    fn preimage_to_fs(
        &self,
        params: &<<Self::M as PolyMatrix>::P as crate::poly::Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
        preimage_id: &str,
    ) -> Vec<std::path::PathBuf> {
        let d = public_matrix.row_size();
        let n = params.ring_dimension() as usize;
        let k = params.modulus_bits();
        let s = SPECTRAL_CONSTANT *
            3.0 *
            SIGMA *
            SIGMA *
            (((d * n * k) as f64).sqrt() + ((2 * n) as f64).sqrt() + 4.7);
        todo!()
    }

    fn preimage_from_fs(
        params: &<<Self::M as PolyMatrix>::P as crate::poly::Poly>::Params,
        preimages_paths: &[PathBuf],
    ) -> Self::M {
        todo!()
    }
}

// A function corresponding to lines 260-266 in trapdoor-dcrtpoly.cpp and the `GaussSampGqArbBase`
// function provided by OpenFHE.
fn gauss_samp_gq_arb_base(
    params: &DCRTPolyParams,
    tower_idx: usize,
    syndrome: &DCRTPoly,
    c: f64,
    dgg: f64,
) -> I64Matrix {
    todo!()
}
