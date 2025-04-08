use crate::{
    parallel_iter,
    poly::{
        dcrt::{
            matrix::disk::{i64_matrix::I64MatrixParams, I64Matrix},
            DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, FinRingElem,
        },
        Poly, PolyParams,
    },
};
use openfhe::ffi::GenerateIntegerKarney;
use rand::{rng, Rng};
use rand_distr::Uniform;
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
    let mut vec = I64Matrix::new_empty(&I64MatrixParams, size, 1);
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
                        val = find_in_vec(m_table, tmp) as i64 * sign;
                    }
                    vec![val]
                })
                .collect::<Vec<Vec<i64>>>()
        };
        vec.replace_entries(0..size, 0..1, f);
    }
    vec
}

pub(crate) fn split_int64_mat_to_elems(
    matrix: &I64Matrix,
    params: &DCRTPolyParams,
) -> DCRTPolyMatrix {
    let n = params.ring_dimension() as usize;
    let nrow = matrix.nrow / n;
    let ncol = matrix.ncol;
    let mut poly_vec = DCRTPolyMatrix::new_empty(params, nrow, ncol);
    let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
        let col_offsets_len = col_offsets.len();
        let i64_values =
            &matrix.block_entries(row_offsets.start * n..row_offsets.end * n, col_offsets);
        parallel_iter!(0..row_offsets.len())
            .map(|i| {
                parallel_iter!(0..col_offsets_len)
                    .map(|j| {
                        let coeffs = i64_values[i * n..(i + 1) * n]
                            .par_iter()
                            .map(|vec| FinRingElem::from_int64(vec[j], params.modulus()))
                            .collect::<Vec<_>>();
                        DCRTPoly::from_coeffs(params, &coeffs)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<Vec<DCRTPoly>>>()
    };
    poly_vec.replace_entries(0..nrow, 0..1, f);
    poly_vec
}

pub(crate) fn split_int64_mat_alt_to_elems(
    matrix: &Vec<Vec<i64>>,
    params: &DCRTPolyParams,
) -> Vec<Vec<DCRTPoly>> {
    let nrow = matrix.len();
    parallel_iter!(0..nrow)
        .map(|i| {
            let coeffs = matrix[i]
                .par_iter()
                .map(|x| FinRingElem::from_int64(*x, params.modulus()))
                .collect::<Vec<_>>();
            vec![DCRTPoly::from_coeffs(params, &coeffs)]
        })
        .collect::<Vec<Vec<DCRTPoly>>>()
}
