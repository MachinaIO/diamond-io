use crate::{
    parallel_iter,
    poly::{
        dcrt::{
            matrix::{
                block_size, dcrt_poly_matrix::CppMatrix, i64_matrix::I64MatrixParams, I64Matrix,
            },
            sampler::DCRTPolyUniformSampler,
            DCRTPolyMatrix, DCRTPolyParams,
        },
        sampler::{DistType, PolyUniformSampler},
        PolyMatrix, PolyParams,
    },
    utils::debug_mem,
};
use openfhe::ffi::{ExtractMatrixCols, FormatMatrixCoefficient, SampleP1ForPertMat};
use rayon::iter::ParallelIterator;
pub use sampler::DCRTPolyTrapdoorSampler;
use std::{cmp::min, ops::Range, sync::Arc};
use utils::{gen_dgg_int_vec, gen_int_karney, split_int64_mat_to_elems};

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
        let log_base_q = params.modulus_digits();
        let dist = DistType::GaussDist { sigma };
        let r = uniform_sampler.sample_uniform(params, size, size * log_base_q, dist);
        let e = uniform_sampler.sample_uniform(params, size, size * log_base_q, dist);
        Self { r, e }
    }

    pub fn sample_pert_square_mat(
        &self,
        s: f64,
        c: f64,
        dgg: f64,
        dgg_large_params: (Option<f64>, f64, Option<&[f64]>),
        peikert: bool,
        total_ncol: usize,
    ) -> DCRTPolyMatrix {
        let r = &self.r;
        let e = &self.e;
        let params = &r.params;
        let n = params.ring_dimension() as usize;
        let (d, dk) = r.size();
        let sigma_large = dgg_large_params.1;
        let num_blocks = total_ncol.div_ceil(d);
        let padded_ncol = num_blocks * d;
        let padding_ncol = padded_ncol - total_ncol;
        debug_mem("sample_pert_square_mat parameters computed");
        // for distribution parameters up to the experimentally found threshold, use
        // the Peikert's inversion method otherwise, use Karney's method
        let p2z_vec = if sigma_large > KARNEY_THRESHOLD {
            let mut matrix = I64Matrix::new_empty(&I64MatrixParams, n * dk, padded_ncol);
            let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<i64>> {
                parallel_iter!(row_offsets)
                    .map(|_| {
                        parallel_iter!(col_offsets.clone())
                            .map(|_| gen_int_karney(0.0, sigma_large))
                            .collect()
                    })
                    .collect()
            };
            matrix.replace_entries(0..n * dk, 0..padded_ncol, f);
            matrix
        } else {
            let dgg_vectors = gen_dgg_int_vec(
                n * dk * padded_ncol,
                peikert,
                dgg_large_params.0.unwrap(),
                dgg_large_params.1,
                dgg_large_params.2.unwrap(),
            );
            let vecs = parallel_iter!(0..n * dk)
                .map(|i| {
                    dgg_vectors.slice(i * padded_ncol, (i + 1) * padded_ncol, 0, 1).transpose()
                })
                .collect::<Vec<_>>();
            vecs[0].concat_rows(&vecs[1..].iter().collect::<Vec<_>>())
        };
        debug_mem("p2z_vec generated");
        // create a matrix of d*k x padded_ncol ring elements in coefficient representation
        let p2 = split_int64_mat_to_elems(&p2z_vec, params);
        // parallel_iter!(0..padded_ncol)
        //     .map(|i| split_int64_vec_to_elems(&p2z_vec.slice(0, n * dk, i, i + 1), params))
        //     .collect::<Vec<_>>();
        // debug_mem("p2_vecs generated");
        // let p2 = p2_vecs[0].concat_columns(&p2_vecs[1..].iter().collect::<Vec<_>>());
        debug_mem("p2 generated");
        let a_mat = r.clone() * r.transpose(); // d x d
        let b_mat = r.clone() * e.transpose(); // d x d
        let d_mat = e.clone() * e.transpose(); // d x d
        debug_mem("a_mat, b_mat, d_mat generated");
        let re = r.concat_rows(&[e]);
        debug_mem("re generated");
        let tp2 = re * &p2;
        debug_mem("tp2 generated");
        let p1 = sample_p1_for_pert_mat(a_mat, b_mat, d_mat, tp2, params, c, s, dgg, padded_ncol);
        debug_mem("p1 generated");
        let mut p = p1.concat_rows(&[&p2]);
        debug_mem("p1 and p2 concatenated");
        if padding_ncol > 0 {
            p = p.slice_columns(0, total_ncol);
        }
        p
    }
}

// A function corresponding to lines 425-473 (except for the line 448) in the `SamplePertSquareMat` function in the trapdoor.h of OpenFHE. https://github.com/openfheorg/openfhe-development/blob/main/src/core/include/lattice/trapdoor.h#L425-L473
fn sample_p1_for_pert_mat(
    a_mat: DCRTPolyMatrix,
    b_mat: DCRTPolyMatrix,
    d_mat: DCRTPolyMatrix,
    tp2: DCRTPolyMatrix,
    params: &DCRTPolyParams,
    c: f64,
    s: f64,
    dgg_stddev: f64,
    padded_ncol: usize,
) -> DCRTPolyMatrix {
    let n = params.ring_dimension();
    let depth = params.crt_depth();
    let k_res = params.crt_bits();
    let block_size = block_size();
    let num_blocks = padded_ncol.div_ceil(block_size);
    let num_threads = rayon::current_num_threads();
    let num_threads_for_cpp = num_threads.div_ceil(num_blocks);
    debug_mem("sample_p1_for_pert_square_mat parameters computed");
    let mut a_mat = a_mat.to_cpp_matrix_ptr();
    FormatMatrixCoefficient(a_mat.inner.as_mut().unwrap());
    let a_mat_arc = Arc::new(a_mat);
    let mut b_mat = b_mat.to_cpp_matrix_ptr();
    FormatMatrixCoefficient(b_mat.inner.as_mut().unwrap());
    let b_mat_arc = Arc::new(b_mat);
    let mut d_mat = d_mat.to_cpp_matrix_ptr();
    FormatMatrixCoefficient(d_mat.inner.as_mut().unwrap());
    let d_mat_arc = Arc::new(d_mat);
    debug_mem("a_mat, b_mat, d_mat are converted to cpp matrices");

    let p1_mat_blocks = parallel_iter!(0..num_blocks)
        .map(|i| {
            let end_col = min((i + 1) * block_size, padded_ncol);
            let mut tp2 = tp2.slice_columns(i * block_size, end_col).to_cpp_matrix_ptr();
            FormatMatrixCoefficient(tp2.inner.as_mut().unwrap());
            let tp2_arc = Arc::new(tp2);
            debug_mem("tp2 is converted to cpp matrices");
            let ncol = end_col - i * block_size;
            let ncol_per_thread = ncol.div_ceil(num_threads_for_cpp);
            let p1_blocks = parallel_iter!(0..ncol.div_ceil(ncol_per_thread))
                .map(|j| {
                    let start_col = j * ncol_per_thread;
                    let end_col = min((j + 1) * ncol_per_thread, ncol);
                    let tp2_cols =
                        ExtractMatrixCols(&Arc::clone(&tp2_arc).as_ref().inner, start_col, end_col);
                    debug_mem("extracting rows from tp2");
                    let cpp_matrix = SampleP1ForPertMat(
                        &Arc::clone(&a_mat_arc).as_ref().inner,
                        &Arc::clone(&b_mat_arc).as_ref().inner,
                        &Arc::clone(&d_mat_arc).as_ref().inner,
                        &tp2_cols,
                        n,
                        depth,
                        k_res,
                        end_col - start_col,
                        c,
                        s,
                        dgg_stddev,
                    );
                    debug_mem("SampleP1ForPertSquareMat called");
                    DCRTPolyMatrix::from_cpp_matrix_ptr(params, &CppMatrix::new(cpp_matrix))
                })
                .collect::<Vec<_>>();
            p1_blocks[0].concat_columns(&p1_blocks[1..].iter().collect::<Vec<_>>())
        })
        .collect::<Vec<_>>();

    p1_mat_blocks[0].concat_columns(&p1_mat_blocks[1..].iter().collect::<Vec<_>>())
}
