use crate::{
    parallel_iter,
    poly::{
        dcrt::{
            matrix::{i64_matrix::I64MatrixParams, I64Matrix},
            DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, FinRingElem,
        },
        Poly, PolyParams,
    },
};
use openfhe::ffi::{DCRTPolyGadgetVector, GenerateIntegerKarney};
use rand::{rng, Rng};
use rand_distr::Uniform;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Range;

pub(crate) fn gen_int_karney(mean: f64, stddev: f64) -> i64 {
    GenerateIntegerKarney(mean, stddev)
}

fn find_in_vec(vec: &[f64], search: f64) -> u32 {
    // binary search to find the position of a value
    let pos = vec.partition_point(|&x| x < search);
    if pos < vec.len() {
        // returns 1-indexed position
        (pos + 1) as u32
    } else {
        panic!("Value not found: {}", search)
    }
}

pub(crate) fn gen_dgg_int_vec(
    size: usize,
    peikert: bool,
    m_a: f64,
    m_std: f64,
    m_table: &[f64],
) -> I64Matrix {
    let mut vec = I64Matrix::zero(&I64MatrixParams, size, 1);
    if !peikert {
        // Use Karney's method
        let f = |row_offsets: Range<usize>, _: Range<usize>| -> Vec<Vec<i64>> {
            parallel_iter!(row_offsets)
                .map(|_| vec![gen_int_karney(0.0f64, m_std)])
                .collect::<Vec<Vec<i64>>>()
        };
        vec.replace_entries(0..size, 0..1, f);
    } else {
        // Use Peikert's algorithm
        let distribution = Uniform::new(0.0f64, 1.0f64).unwrap();
        let f = |row_offsets: Range<usize>, _: Range<usize>| -> Vec<Vec<i64>> {
            parallel_iter!(row_offsets)
                .map(|_| {
                    let mut rng = rng();
                    let seed: f64 = rng.sample(distribution) - 0.5f64;
                    let tmp = seed.abs() - m_a / 2.0f64;
                    let mut val = 0;

                    if tmp > 0.0f64 {
                        let sign = if seed > 0.0f64 { 1 } else { -1 };
                        val = find_in_vec(&m_table, tmp) as i64 * sign;
                    }
                    vec![val]
                })
                .collect::<Vec<Vec<i64>>>()
        };
        vec.replace_entries(0..size, 0..1, f);
    }
    vec
}

pub(crate) fn split_int64_vec_to_elems(vec: &I64Matrix, params: &DCRTPolyParams) -> DCRTPolyMatrix {
    debug_assert_eq!(vec.ncol, 1, "Matrix must be a column vector");
    let n = params.ring_dimension() as usize;
    let nrow = vec.nrow / n;
    let mut poly_vec = DCRTPolyMatrix::zero(params, nrow, 1);
    let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
        debug_assert_eq!(col_offsets.len(), 1, "Matrix must be a column vector");
        let i64_values = &vec
            .block_entries(row_offsets.start * n..row_offsets.end * n, col_offsets)
            .into_iter()
            .map(|vec| vec[0])
            .collect::<Vec<_>>();
        parallel_iter!(0..row_offsets.len())
            .map(|i| {
                let coeffs = i64_values[i * n..(i + 1) * n]
                    .iter()
                    .map(|x| FinRingElem::from_int64(*x, params.modulus()))
                    .collect::<Vec<_>>();
                vec![DCRTPoly::from_coeffs(params, &coeffs)]
            })
            .collect::<Vec<Vec<DCRTPoly>>>()
    };
    poly_vec.replace_entries(0..nrow, 0..1, f);
    poly_vec
}

pub(crate) fn split_int64_vec_alt_to_elems(
    vec: &I64Matrix,
    params: &DCRTPolyParams,
) -> DCRTPolyMatrix {
    let n = params.ring_dimension() as usize;
    debug_assert_eq!(vec.ncol, n, "Matrix must have n columns");
    let nrow = vec.nrow;
    let mut poly_vec = DCRTPolyMatrix::zero(params, nrow, 1);
    let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
        debug_assert_eq!(col_offsets.len(), 1, "Matrix must be a column vector");
        let i64_values = &vec.block_entries(row_offsets.clone(), 0..n);
        parallel_iter!(0..row_offsets.len())
            .map(|i| {
                let coeffs = i64_values[i]
                    .iter()
                    .map(|x| FinRingElem::from_int64(*x, params.modulus()))
                    .collect::<Vec<_>>();
                vec![DCRTPoly::from_coeffs(params, &coeffs)]
            })
            .collect::<Vec<Vec<DCRTPoly>>>()
    };
    poly_vec.replace_entries(0..nrow, 0..1, f);
    poly_vec
}

pub(crate) fn gen_dcrt_gadget_vector(params: &DCRTPolyParams) -> DCRTPolyMatrix {
    let g_vec_cpp = DCRTPolyGadgetVector(
        params.ring_dimension(),
        params.crt_depth(),
        params.crt_bits(),
        params.modulus_bits(),
        2,
    );
    DCRTPolyMatrix::from_cpp_matrix_ptr(params, g_vec_cpp)
}
