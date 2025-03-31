use super::{block_offsets, MmapMatrix, MmapMatrixElem, MmapMatrixParams};
use crate::parallel_iter;
use crate::poly::dcrt::FinRingElem;
use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams},
    Poly, PolyMatrix, PolyParams,
};
use itertools::Itertools;
use num_bigint::BigInt;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Range;

impl MmapMatrixParams for DCRTPolyParams {
    fn entry_size(&self) -> usize {
        let log_q_bytes = self.modulus_bits().div_ceil(8);
        let dim = self.ring_dimension() as usize;
        dim * log_q_bytes
    }
}

impl MmapMatrixElem for DCRTPoly {
    type Params = DCRTPolyParams;

    fn zero(params: &Self::Params) -> Self {
        <Self as Poly>::const_zero(params)
    }
    fn one(params: &Self::Params) -> Self {
        <Self as Poly>::const_one(params)
    }
    fn from_bytes_to_elem(params: &Self::Params, bytes: &[u8]) -> Self {
        <Self as Poly>::from_bytes(params, bytes)
    }

    fn from_elem_to_bytes(&self) -> Vec<u8> {
        self.to_bytes()
    }
}

pub type DCRTPolyMatrix = MmapMatrix<DCRTPoly>;

impl PolyMatrix for DCRTPolyMatrix {
    type P = DCRTPoly;

    fn from_poly_vec(params: &DCRTPolyParams, vec: Vec<Vec<DCRTPoly>>) -> Self {
        let nrow = vec.len();
        let ncol = vec[0].len();
        let mut matrix = Self::new_empty(params, nrow, ncol);
        let vec = &vec;
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<Self::P>> {
            let mut new_entries =
                vec![vec![DCRTPoly::const_zero(params); col_offsets.len()]; row_offsets.len()];
            let row_start = row_offsets.start;
            for i in row_offsets {
                new_entries[i - row_start].clone_from_slice(&vec[i][col_offsets.clone()]);
            }
            new_entries
        };
        matrix.replace_entries(0..nrow, 0..ncol, f);
        matrix
    }

    fn entry(&self, i: usize, j: usize) -> Self::P {
        self.entry(i, j)
    }

    fn get_row(&self, i: usize) -> Vec<Self::P> {
        self.get_row(i)
    }

    fn get_column(&self, j: usize) -> Vec<Self::P> {
        self.get_column(j)
    }

    fn size(&self) -> (usize, usize) {
        self.size()
    }

    fn slice(&self, row_start: usize, row_end: usize, col_start: usize, col_end: usize) -> Self {
        self.slice(row_start, row_end, col_start, col_end)
    }

    fn zero(params: &<Self::P as Poly>::Params, nrow: usize, ncol: usize) -> Self {
        Self::zero(params, nrow, ncol)
    }

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        Self::identity(params, size, scalar)
    }

    fn transpose(&self) -> Self {
        self.transpose()
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[&Self]) -> Self {
        self.concat_columns(others)
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[&Self]) -> Self {
        self.concat_rows(others)
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[&Self]) -> Self {
        self.concat_diag(others)
    }

    fn tensor(&self, other: &Self) -> Self {
        self.tensor(other)
    }

    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self {
        let bit_length = params.modulus_bits();
        let mut poly_vec = Vec::with_capacity(bit_length);
        for i in 0u32..(bit_length as u32) {
            let value = BigInt::from(2).pow(i);
            poly_vec.push(DCRTPoly::from_const(params, &FinRingElem::new(value, params.modulus())));
        }
        let gadget_vector = Self::from_poly_vec(params, vec![poly_vec]);
        let identity = DCRTPolyMatrix::identity(params, size, None);
        identity.tensor(&gadget_vector)
    }

    fn decompose(&self) -> Self {
        let bit_length = self.params.modulus_bits();
        let new_nrow = self.nrow * bit_length;
        let new_matrix = Self::new_empty(&self.params, new_nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let block_row_len = next_block_row_idx - cur_block_row_idx;
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        let mut new_entries =
                            vec![
                                vec![DCRTPoly::const_zero(&self.params); block_col_len];
                                block_row_len * bit_length
                            ];
                        for j in 0..block_col_len {
                            for i in 0..block_row_len {
                                let poly = &self_block_polys[i][j];
                                let decomposed_polys = poly.decompose(&self.params);
                                for (k, poly) in decomposed_polys.into_iter().enumerate() {
                                    new_entries[i * bit_length + k][j] = poly;
                                }
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                cur_block_row_idx * bit_length..next_block_row_idx * bit_length,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_entries,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }

    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self {
        let mut new_matrix = Self::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<Self::P>> {
            let self_block_polys = self.block_entries(row_offsets, col_offsets);
            self_block_polys
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|poly| poly.modulus_switch(&self.params, new_modulus.clone()))
                        .collect_vec()
                })
                .collect_vec()
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }

    fn mul_tensor_identity(&self, other: &Self, identity_size: usize) -> Self {
        debug_assert_eq!(self.ncol, other.nrow * identity_size);

        let slice_width = other.nrow;
        let mut slice_results = Vec::with_capacity(identity_size);

        for i in 0..identity_size {
            let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
            slice_results.push(slice * other);
        }
        slice_results[0].clone().concat_columns(&slice_results[1..].iter().collect::<Vec<_>>())
    }

    fn mul_tensor_identity_decompose(&self, other: &Self, identity_size: usize) -> Self {
        let log_q = self.params.modulus_bits();
        debug_assert_eq!(self.ncol, other.nrow * identity_size * log_q);
        let slice_width = other.nrow * log_q;
        let mut output = vec![DCRTPolyMatrix::zero(&self.params, 0, 0); other.ncol * identity_size];

        for j in 0..other.ncol {
            let jth_col_m_decompose = other.get_column_matrix_decompose(j);
            for i in 0..identity_size {
                let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
                output[i * other.ncol + j] = slice * &jth_col_m_decompose;
            }
        }

        output[0].clone().concat_columns(&output[1..].iter().collect::<Vec<_>>())
    }

    fn get_column_matrix_decompose(&self, j: usize) -> Self {
        Self::from_poly_vec(
            &self.params,
            self.get_column(j).into_iter().map(|poly| vec![poly]).collect(),
        )
        .decompose()
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use std::sync::Arc;

    use num_bigint::BigUint;

    use super::*;
    use crate::poly::{
        dcrt::{DCRTPolyParams, DCRTPolyUniformSampler},
        sampler::PolyUniformSampler,
    };

    #[test]
    fn test_matrix_gadget_matrix() {
        let params = DCRTPolyParams::default();
        let size = 3;
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, size);
        assert_eq!(gadget_matrix.size().0, size);
        assert_eq!(gadget_matrix.size().1, size * params.modulus_bits());
    }

    #[test]
    fn test_matrix_decompose() {
        let params = DCRTPolyParams::default();
        let bit_length = params.modulus_bits();

        // Create a simple 2x8 matrix with some non-zero values
        let mut matrix_vec = Vec::with_capacity(2);
        let value = FinRingElem::new(5u32, params.modulus());

        // Create first row
        let mut row1 = Vec::with_capacity(8);
        row1.push(DCRTPoly::from_const(&params, &value));
        for _ in 1..8 {
            row1.push(DCRTPoly::const_zero(&params));
        }

        // Create second row
        let mut row2 = Vec::with_capacity(8);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::from_const(&params, &value));
        for _ in 2..8 {
            row2.push(DCRTPoly::const_zero(&params));
        }

        matrix_vec.push(row1);
        matrix_vec.push(row2);

        let matrix = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
        assert_eq!(matrix.size().0, 2);
        assert_eq!(matrix.size().1, 8);

        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, 2);
        assert_eq!(gadget_matrix.size().0, 2);
        assert_eq!(gadget_matrix.size().1, 2 * bit_length);

        let decomposed = matrix.decompose();
        assert_eq!(decomposed.size().0, 2 * bit_length);
        assert_eq!(decomposed.size().1, 8);

        let expected_matrix = gadget_matrix * decomposed;
        assert_eq!(expected_matrix.size().0, 2);
        assert_eq!(expected_matrix.size().1, 8);
        assert_eq!(matrix, expected_matrix);
    }

    #[test]
    fn test_matrix_basic_operations() {
        let params = DCRTPolyParams::default();

        // Test zero and identity matrices
        let zero = DCRTPolyMatrix::zero(&params, 2, 2);
        let identity = DCRTPolyMatrix::identity(&params, 2, None);

        // Test matrix creation and equality
        let value = FinRingElem::new(5u32, params.modulus());

        // Create a 2x2 matrix with values at (0,0) and (1,1)
        let matrix_vec = vec![
            vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
            vec![DCRTPoly::const_zero(&params), DCRTPoly::from_const(&params, &value)],
        ];

        let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
        assert_eq!(matrix1.entry(0, 0).coeffs()[0], value);
        let matrix2 = matrix1.clone();
        assert_eq!(matrix1, matrix2);

        // Test addition
        let sum = matrix1.clone() + &matrix2;
        let value_10 = FinRingElem::new(10u32, params.modulus());
        assert_eq!(sum.entry(0, 0).coeffs()[0], value_10);

        // Test subtraction
        let diff = matrix1.clone() - &matrix2;
        assert_eq!(diff, zero);

        // Test multiplication
        let prod = matrix1 * &identity;
        assert_eq!(prod.size(), (2, 2));
        // Check that the product has the same values as the original matrix
        assert_eq!(prod.entry(0, 0).coeffs()[0], value);
        assert_eq!(prod.entry(1, 1).coeffs()[0], value);
    }

    #[test]
    fn test_matrix_concatenation() {
        let params = DCRTPolyParams::default();
        let value = FinRingElem::new(5u32, params.modulus());

        // Create first matrix with value at (0,0)
        let matrix1_vec = vec![
            vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
            vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
        ];

        let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix1_vec);

        // Create second matrix with value at (1,1)
        let matrix2_vec = vec![
            vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
            vec![DCRTPoly::const_zero(&params), DCRTPoly::from_const(&params, &value)],
        ];

        let matrix2 = DCRTPolyMatrix::from_poly_vec(&params, matrix2_vec);

        // Test column concatenation
        let col_concat = matrix1.concat_columns(&[&matrix2]);
        assert_eq!(col_concat.size().0, 2);
        assert_eq!(col_concat.size().1, 4);
        assert_eq!(col_concat.entry(0, 0).coeffs()[0], value);
        assert_eq!(col_concat.entry(1, 3).coeffs()[0], value);

        // Test row concatenation
        let row_concat = matrix1.concat_rows(&[&matrix2]);
        assert_eq!(row_concat.size().0, 4);
        assert_eq!(row_concat.size().1, 2);
        assert_eq!(row_concat.entry(0, 0).coeffs()[0], value);
        assert_eq!(row_concat.entry(3, 1).coeffs()[0], value);

        // Test diagonal concatenation
        let diag_concat = matrix1.concat_diag(&[&matrix2]);
        assert_eq!(diag_concat.size().0, 4);
        assert_eq!(diag_concat.size().1, 4);
        assert_eq!(diag_concat.entry(0, 0).coeffs()[0], value);
        assert_eq!(diag_concat.entry(3, 3).coeffs()[0], value);
    }

    #[test]
    fn test_matrix_tensor_product() {
        let params = DCRTPolyParams::default();
        let value = FinRingElem::new(5u32, params.modulus());

        // Create first matrix with value at (0,0)
        let matrix1_vec = vec![
            vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
            vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
        ];

        let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix1_vec);

        // Create second matrix with value at (0,0)
        let matrix2_vec = vec![
            vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
            vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
        ];

        let matrix2 = DCRTPolyMatrix::from_poly_vec(&params, matrix2_vec);

        let tensor = matrix1.tensor(&matrix2);
        assert_eq!(tensor.size().0, 4);
        assert_eq!(tensor.size().1, 4);

        // Check that the (0,0) element is the product of the (0,0) elements
        let value_25 = FinRingElem::new(25u32, params.modulus());
        assert_eq!(tensor.entry(0, 0).coeffs()[0], value_25);
    }

    #[test]
    fn test_matrix_modulus_switch() {
        let params = DCRTPolyParams::default();

        let value00 = FinRingElem::new(1023782870921908217643761278891282178u128, params.modulus());
        let value01 = FinRingElem::new(8179012198875468938912873783289218738u128, params.modulus());
        let value10 = FinRingElem::new(2034903202902173762872163465127672178u128, params.modulus());
        let value11 = FinRingElem::new(1990091289902891278121564387120912660u128, params.modulus());

        let matrix_vec = vec![
            vec![DCRTPoly::from_const(&params, &value00), DCRTPoly::from_const(&params, &value01)],
            vec![DCRTPoly::from_const(&params, &value10), DCRTPoly::from_const(&params, &value11)],
        ];

        let matrix = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
        let new_modulus = Arc::new(BigUint::from(2u32));
        let switched = matrix.modulus_switch(&new_modulus);

        // Although the value becomes less than the new modulus, the set modulus is still the same
        assert_eq!(switched.params.modulus(), params.modulus());

        let new_value00 = value00.modulus_switch(new_modulus.clone());
        let new_value01 = value01.modulus_switch(new_modulus.clone());
        let new_value10 = value10.modulus_switch(new_modulus.clone());
        let new_value11 = value11.modulus_switch(new_modulus.clone());

        let expected_vec = vec![
            vec![
                DCRTPoly::from_const(&params, &new_value00),
                DCRTPoly::from_const(&params, &new_value01),
            ],
            vec![
                DCRTPoly::from_const(&params, &new_value10),
                DCRTPoly::from_const(&params, &new_value11),
            ],
        ];

        let expected = DCRTPolyMatrix::from_poly_vec(&params, expected_vec);
        assert_eq!(switched, expected);
    }

    #[test]
    #[should_panic(expected = "Addition requires matrices of same dimensions")]
    #[cfg(debug_assertions)]
    fn test_matrix_addition_mismatch() {
        let params = DCRTPolyParams::default();
        let matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        let matrix2 = DCRTPolyMatrix::zero(&params, 2, 3);
        let _sum = matrix1 + matrix2;
    }

    #[test]
    #[should_panic(expected = "Multiplication condition failed")]
    #[cfg(debug_assertions)]
    fn test_matrix_multiplication_mismatch() {
        let params = DCRTPolyParams::default();
        let matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        let matrix2 = DCRTPolyMatrix::zero(&params, 3, 2);
        let _prod = matrix1 * matrix2;
    }

    #[test]
    fn test_matrix_mul_tensor_identity_simple() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix S (2x20)
        let s = sampler.sample_uniform(&params, 2, 20, crate::poly::sampler::DistType::FinRingDist);
        // Create 'other' matrix (5x7)
        let other =
            sampler.sample_uniform(&params, 5, 7, crate::poly::sampler::DistType::FinRingDist);
        // Perform S * (I_4 ⊗ other)
        let result = s.mul_tensor_identity(&other, 4);

        // Check dimensions
        assert_eq!(result.size().0, 2);
        assert_eq!(result.size().1, 28);

        let identity = DCRTPolyMatrix::identity(&params, 4, None);
        // Check result
        let expected_result = s * (identity.tensor(&other));

        assert_eq!(expected_result.size().0, 2);
        assert_eq!(expected_result.size().1, 28);
        assert_eq!(result, expected_result)
    }

    #[test]
    fn test_matrix_mul_tensor_identity_decompose_naive() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix S (2x2516)
        let s =
            sampler.sample_uniform(&params, 2, 2516, crate::poly::sampler::DistType::FinRingDist);

        // Create 'other' matrix (2x13)
        let other =
            sampler.sample_uniform(&params, 2, 13, crate::poly::sampler::DistType::FinRingDist);

        // Decompose 'other' matrix
        let other_decompose = other.decompose();
        // Perform S * (I_37 ⊗ G^-1(other))
        let result: DCRTPolyMatrix = s.mul_tensor_identity(&other_decompose, 37);
        // Check dimensions
        assert_eq!(result.size().0, 2);
        assert_eq!(result.size().1, 481);

        // Check result
        let tensor = identity_tensor_matrix(37, &other_decompose);
        let expected_result = s * tensor;

        assert_eq!(expected_result.size().0, 2);
        assert_eq!(expected_result.size().1, 481);
        assert_eq!(result, expected_result)
    }

    #[test]
    fn test_matrix_mul_tensor_identity_decompose_optimal() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix S (2x2516)
        let s =
            sampler.sample_uniform(&params, 2, 2516, crate::poly::sampler::DistType::FinRingDist);

        // Create 'other' matrix (2x13)
        let other =
            sampler.sample_uniform(&params, 2, 13, crate::poly::sampler::DistType::FinRingDist);

        // Perform S * (I_37 ⊗ G^-1(other))
        let result: DCRTPolyMatrix = s.mul_tensor_identity_decompose(&other, 37);

        // Check dimensions
        assert_eq!(result.size().0, 2);
        assert_eq!(result.size().1, 481);

        // Check result
        let decomposed = other.decompose();
        let tensor = identity_tensor_matrix(37, &decomposed);
        let expected_result_1 = s.clone() * tensor;
        let expected_result_2 = s.mul_tensor_identity(&decomposed, 37);
        assert_eq!(expected_result_1, expected_result_2);

        assert_eq!(expected_result_1.size().0, 2);
        assert_eq!(expected_result_1.size().1, 481);

        assert_eq!(expected_result_2.size().0, 2);
        assert_eq!(expected_result_2.size().1, 481);

        assert_eq!(result, expected_result_1);
        assert_eq!(result, expected_result_2);
    }

    fn identity_tensor_matrix(identity_size: usize, matrix: &DCRTPolyMatrix) -> DCRTPolyMatrix {
        let mut others = vec![];
        for _ in 1..identity_size {
            others.push(matrix);
        }
        matrix.concat_diag(&others[..])
    }
}
