use crate::{
    parallel_iter,
    poly::{
        dcrt::{
            matrix::{i64_matrix::I64MatrixParams, I64Matrix},
            sampler::DCRTPolyUniformSampler,
            DCRTPolyMatrix, DCRTPolyParams,
        },
        sampler::{DistType, PolyUniformSampler},
        PolyParams,
    },
    utils::debug_mem,
};
use openfhe::ffi::SampleP1ForPertSquareMat;
use rayon::iter::ParallelIterator;
pub use sampler::DCRTPolyTrapdoorSampler;
use std::ops::Range;
use utils::{gen_dgg_int_vec, gen_int_karney, split_int64_vec_to_elems};

pub mod sampler;
pub mod utils;

pub(crate) const KARNEY_THRESHOLD: f64 = 300.0;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DCRTTrapdoor {
    pub r: DCRTPolyMatrix,
    pub e: DCRTPolyMatrix,
}

impl DCRTTrapdoor {
    pub fn new(params: &DCRTPolyParams, size: usize, sigma: f64) -> Self {
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let log_q = params.modulus_bits();
        let dist = DistType::GaussDist { sigma };
        let r = uniform_sampler.sample_uniform(params, size, size * log_q, dist);
        let e = uniform_sampler.sample_uniform(params, size, size * log_q, dist);
        Self { r, e }
    }

    pub fn sample_pert_square_mat(
        &self,
        s: f64,
        c: f64,
        dgg: f64,
        dgg_large_params: (f64, f64, &[f64]),
        peikert: bool,
    ) -> DCRTPolyMatrix {
        let r = &self.r;
        let e = &self.e;
        let params = &r.params;
        let n = params.ring_dimension() as usize;
        let (d, dk) = r.size();
        let sigma_large = (s * s - c * c).sqrt();
        debug_mem("sample_pert_square_mat parameters computed");
        // for distribution parameters up to the experimentally found threshold, use
        // the Peikert's inversion method otherwise, use Karney's method
        let p2z_vec = if sigma_large > KARNEY_THRESHOLD {
            let mut matrix = I64Matrix::new_empty(&I64MatrixParams, n * dk, d);
            let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<i64>> {
                parallel_iter!(row_offsets)
                    .map(|_| {
                        parallel_iter!(col_offsets.clone())
                            .map(|_| gen_int_karney(0.0, sigma_large))
                            .collect()
                    })
                    .collect()
            };
            matrix.replace_entries(0..n * dk, 0..d, f);
            matrix
        } else {
            let dgg_vectors = gen_dgg_int_vec(
                d,
                peikert,
                dgg_large_params.0,
                dgg_large_params.1,
                dgg_large_params.2,
            );
            let vecs = parallel_iter!(0..n * dk)
                .map(|i| dgg_vectors.slice(i * d, (i + 1) * d, 0, 1))
                .collect::<Vec<_>>();
            vecs[0].concat_columns(&vecs[1..].iter().collect::<Vec<_>>())
        };
        debug_mem("p2z_vec generated");
        // create a matrix of d*k x d ring elements in coefficient representation
        let p2_vecs = parallel_iter!(0..d)
            .map(|i| split_int64_vec_to_elems(&p2z_vec.slice(0, n * dk, i, i + 1), params))
            .collect::<Vec<_>>();
        let p2 = p2_vecs[0].concat_columns(&p2_vecs[1..].iter().collect::<Vec<_>>());
        debug_mem("p2 generated");
        let a_mat = r.clone() * r.transpose(); // d * d
        let b_mat = r.clone() * e.transpose(); // d * d
        let d_mat = e.clone() * e.transpose(); // d * d
        let tp2 = r.concat_rows(&[e]) * &p2;
        let p1 = sample_p1_for_pert_square_mat(a_mat, b_mat, d_mat, tp2, params, c, s, dgg);
        debug_mem("p1 generated");
        p1.concat_rows(&[&p2])
    }
}

// A function to generate `p1` for the `sample_pert_square_mat` function, corresponding to lines 425-473 (except for the line 448) in the `SamplePertSquareMat` function in the trapdoor.h of OpenFHE. https://github.com/openfheorg/openfhe-development/blob/main/src/core/include/lattice/trapdoor.h#L425-L473
fn sample_p1_for_pert_square_mat(
    a_mat: DCRTPolyMatrix,
    b_mat: DCRTPolyMatrix,
    d_mat: DCRTPolyMatrix,
    tp2: DCRTPolyMatrix,
    params: &DCRTPolyParams,
    c: f64,
    s: f64,
    dgg_stddev: f64,
) -> DCRTPolyMatrix {
    let n = params.ring_dimension();
    let depth = params.crt_depth();
    let k_res = params.modulus_bits() / depth;
    let mut a = a_mat.to_cpp_matrix_ptr();
    let mut b = b_mat.to_cpp_matrix_ptr();
    let mut d = d_mat.to_cpp_matrix_ptr();
    let mut tp2 = tp2.to_cpp_matrix_ptr();
    let p1_mat = SampleP1ForPertSquareMat(
        a.as_mut().unwrap(),
        b.as_mut().unwrap(),
        d.as_mut().unwrap(),
        tp2.as_mut().unwrap(),
        n,
        depth,
        k_res,
        c,
        s,
        dgg_stddev,
    );

    DCRTPolyMatrix::from_cpp_matrix_ptr(params, p1_mat)
}
