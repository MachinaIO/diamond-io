use crate::{
    parallel_iter,
    poly::{
        dcrt::{cpp_matrix::CppMatrix, DCRTPoly, DCRTPolyParams},
        MatrixElem, MatrixParams, Poly, PolyMatrix, PolyParams,
    },
    utils::{block_size, debug_mem},
};
use itertools::Itertools;
use openfhe::ffi::{DCRTPolyGadgetVector, MatrixGen, SetMatrixElement};
use rayon::prelude::*;
use std::{ops::Range, path::Path};

use super::base::BaseMatrix;

#[cfg(feature = "disk")]
use super::base::disk::{block_offsets, map_file_mut};

impl MatrixParams for DCRTPolyParams {
    fn entry_size(&self) -> usize {
        let log_q_bytes = self.modulus_bits().div_ceil(8);
        let dim = self.ring_dimension() as usize;
        dim * log_q_bytes
    }
}

impl MatrixElem for DCRTPoly {
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

    fn as_elem_to_bytes(&self) -> Vec<u8> {
        self.to_bytes()
    }
}

pub type DCRTPolyMatrix = BaseMatrix<DCRTPoly>;

impl PolyMatrix for DCRTPolyMatrix {
    type P = DCRTPoly;

    fn from_poly_vec(params: &DCRTPolyParams, vec: Vec<Vec<DCRTPoly>>) -> Self {
        let nrow = vec.len();
        let ncol = vec[0].len();
        let mut matrix = Self::new_empty(params, nrow, ncol);
        let vec = &vec;
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<Self::P>> {
            row_offsets.into_iter().map(|i| vec[i][col_offsets.clone()].to_vec()).collect()
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
        let gadget_vector = Self::gadget_vector(params);
        debug_assert_eq!(gadget_vector.col_size(), params.modulus_digits());
        gadget_vector.concat_diag(&vec![&gadget_vector; size - 1])
    }

    fn decompose(&self) -> Self {
        let base_bits = self.params.base_bits();
        let log_base_q =
            self.params.crt_bits().div_ceil(base_bits as usize) * self.params.crt_depth();
        let new_nrow = self.nrow * log_base_q;
        let mut new_matrix = Self::new_empty(&self.params, new_nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
            let nrow = row_offsets.len();
            let new_nrow = row_offsets.len() * log_base_q;
            let ncol = col_offsets.len();
            let entries = self.block_entries(row_offsets, col_offsets);
            let decomposed_entries: Vec<Vec<Vec<DCRTPoly>>> = parallel_iter!(0..nrow)
                .map(|i| {
                    parallel_iter!(0..ncol)
                        .map(|j| self.dcrt_decompose_poly(&entries[i][j], base_bits))
                        .collect()
                })
                .collect();
            parallel_iter!(0..new_nrow)
                .map(|idx| {
                    let i = idx / log_base_q;
                    let k = idx % log_base_q;

                    parallel_iter!(0..ncol).map(|j| decomposed_entries[i][j][k].clone()).collect()
                })
                .collect()
        };
        new_matrix.replace_entries_with_expand(0..self.nrow, 0..self.ncol, log_base_q, 1, f);
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

        let slice_results = (0..identity_size)
            .map(|i| {
                let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
                slice * other
            })
            .collect_vec();

        slice_results[0].concat_columns(&slice_results[1..].iter().collect::<Vec<_>>())
    }

    fn mul_tensor_identity_decompose(&self, other: &Self, identity_size: usize) -> Self {
        let log_base_q = self.params.modulus_digits();
        debug_assert_eq!(self.ncol, other.nrow * identity_size * log_base_q);
        let slice_width = other.nrow * log_base_q;

        let output = (0..identity_size)
            .flat_map(|i| {
                debug_mem(format!("mul_tensor_identity_decompose at {i}"));
                let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
                (0..other.ncol).map(move |j| &slice * &other.get_column_matrix_decompose(j))
            })
            .collect_vec();
        debug_mem("mul_tensor_identity_decompose output computed".to_string());
        output[0].concat_columns(&output[1..].iter().collect::<Vec<_>>())
    }

    fn get_column_matrix_decompose(&self, j: usize) -> Self {
        Self::from_poly_vec(
            &self.params,
            self.get_column(j).into_iter().map(|poly| vec![poly]).collect(),
        )
        .decompose()
    }

    #[inline]
    fn read_from_files<P: AsRef<Path> + Send + Sync>(
        params: &<Self::P as Poly>::Params,
        nrow: usize,
        ncol: usize,
        dir_path: P,
        id: &str,
    ) -> Self {
        let block_size = block_size();
        let mut matrix = Self::new_empty(params, nrow, ncol);

        let f = |row_range: Range<usize>, col_range: Range<usize>| -> Vec<Vec<DCRTPoly>> {
            let mut path = dir_path.as_ref().to_path_buf();
            path.push(format!(
                "{}_{}_{}.{}_{}.{}.matrix",
                id, block_size, row_range.start, row_range.end, col_range.start, col_range.end
            ));
            let file = std::fs::File::open(&path)
                .unwrap_or_else(|_| panic!("Failed to read matrix file {path:?}"));
            let mut reader = std::io::BufReader::new(file);
            // Stream decode directly into flat Vec to avoid triple nesting
            let entries_flat: Vec<Vec<u8>> = bincode::decode_from_std_read::<Vec<Vec<u8>>, _, _>(
                &mut reader,
                bincode::config::standard(),
            )
            .unwrap();
            let cols = col_range.len();

            // Keep parallel processing but with flat indexing
            parallel_iter!(0..row_range.len())
                .map(|i| {
                    parallel_iter!(0..cols)
                        .map(|j| {
                            let flat_idx = i * cols + j;
                            DCRTPoly::from_compact_bytes(params, &entries_flat[flat_idx])
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        matrix.replace_entries(0..nrow, 0..ncol, f);
        matrix
    }

    // #[inline]
    // async fn write_to_files<P: AsRef<Path> + Send + Sync>(&self, dir_path: P, id: &str) {
    //     let block_size = block_size();
    //     info!("{block_size}");
    //     #[cfg(feature = "disk")]
    //     let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
    //     #[cfg(not(feature = "disk"))]
    //     let (row_offsets, col_offsets) = (vec![0, self.nrow], vec![0, self.ncol]);
    //     let dir_path = dir_path.as_ref().to_path_buf();

    //     let self_arc = Arc::new(self);
    //     let row_windows = row_offsets.into_iter().tuple_windows().collect_vec();
    //     let futures = row_windows
    //         .into_iter()
    //         .flat_map(|(cur_block_row_idx, next_block_row_idx)| {
    //             let col_windows = col_offsets.clone().into_iter().tuple_windows().collect_vec();
    //             col_windows
    //                 .into_iter()
    //                 .map(|(cur_block_col_idx, next_block_col_idx)| {
    //                     // let id_clone = id.clone();
    //                     let row_range = cur_block_row_idx..next_block_row_idx;
    //                     let col_range = cur_block_col_idx..next_block_col_idx;
    //                     let self_arc = Arc::clone(&self_arc);
    //                     let dir_path = dir_path.clone();
    //                     async move {
    //                         let entries = self_arc
    //                             .as_ref()
    //                             .block_entries(row_range.clone(), col_range.clone());
    //                         let mut path = dir_path;
    //                         path.push(format!(
    //                             "{}_{}_{}.{}_{}.{}.matrix",
    //                             id,
    //                             block_size,
    //                             row_range.start,
    //                             row_range.end,
    //                             col_range.start,
    //                             col_range.end
    //                         ));
    //                         let entries_bytes: Vec<Vec<Vec<u8>>> = entries
    //                             .iter()
    //                             .map(|row| {
    //                                 row.iter().map(|poly| poly.to_compact_bytes()).collect_vec()
    //                             })
    //                             .collect_vec();
    //                         let serialized_data =
    //                             bincode::encode_to_vec(&entries_bytes,
    // bincode::config::standard())                                 .unwrap();
    //                         info!("serialized_len={}", serialized_data.len());
    //                         write(path, &serialized_data).await
    //                     }
    //                 })
    //                 .collect::<Vec<_>>()
    //         })
    //         .collect::<Vec<_>>();
    //     futures::future::try_join_all(futures).await.expect("Failed to write all matrix blocks");
    // }

    fn set_entry(&mut self, i: usize, j: usize, elem: Self::P) {
        #[cfg(not(feature = "disk"))]
        {
            self.inner[i][j] = elem;
        }
        #[cfg(feature = "disk")]
        {
            let entry_size = self.params.entry_size();
            let offset = (i * self.ncol + j) * entry_size;
            let bytes = elem.as_elem_to_bytes();
            unsafe {
                let mut mmap = map_file_mut(&self.file, offset, entry_size);
                mmap.copy_from_slice(&bytes);
            }
        }
    }

    fn to_compact_bytes(&self) -> Vec<u8> {
        let entries = self.block_entries(0..self.nrow, 0..self.ncol);
        let entries_bytes: Vec<Vec<Vec<u8>>> = entries
            .iter()
            .map(|row| row.iter().map(|poly| poly.to_compact_bytes()).collect())
            .collect();

        bincode::encode_to_vec(&entries_bytes, bincode::config::standard())
            .expect("Failed to serialize matrix to compact bytes")
    }

    fn from_compact_bytes(params: &<Self::P as Poly>::Params, bytes: &[u8]) -> Self {
        let entries_bytes: Vec<Vec<Vec<u8>>> =
            bincode::decode_from_slice(bytes, bincode::config::standard())
                .expect("Failed to deserialize matrix from compact bytes")
                .0;

        let nrow = entries_bytes.len();
        let ncol = if nrow > 0 { entries_bytes[0].len() } else { 0 };
        let mut matrix = Self::new_empty(params, nrow, ncol);

        let f = |row_range: Range<usize>, col_range: Range<usize>| -> Vec<Vec<DCRTPoly>> {
            row_range
                .map(|i| {
                    col_range
                        .clone()
                        .map(|j| DCRTPoly::from_compact_bytes(params, &entries_bytes[i][j]))
                        .collect()
                })
                .collect()
        };

        matrix.replace_entries(0..nrow, 0..ncol, f);
        matrix
    }

    fn block_entries(
        &self,
        rows: std::ops::Range<usize>,
        cols: std::ops::Range<usize>,
    ) -> Vec<Vec<Self::P>> {
        // Delegate to the BaseMatrix implementation
        self.block_entries(rows, cols)
    }
}

impl DCRTPolyMatrix {
    pub(crate) fn to_cpp_matrix_ptr(&self) -> CppMatrix {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let mut matrix_ptr = MatrixGen(
            self.params.ring_dimension(),
            self.params.crt_depth(),
            self.params.crt_bits(),
            nrow,
            ncol,
        );
        debug_mem(format!("matrix_ptr MatrixGen row={nrow}, col={ncol}"));
        for i in 0..nrow {
            for j in 0..ncol {
                SetMatrixElement(matrix_ptr.as_mut().unwrap(), i, j, self.entry(i, j).get_poly());
            }
        }
        debug_mem(format!("SetMatrixElement row={nrow}, col={ncol}"));
        CppMatrix::new(&self.params, matrix_ptr)
    }

    pub(crate) fn from_cpp_matrix_ptr(params: &DCRTPolyParams, cpp_matrix: &CppMatrix) -> Self {
        let nrow = cpp_matrix.nrow();
        let ncol = cpp_matrix.ncol();
        let matrix_inner = parallel_iter!(0..nrow)
            .map(|i| parallel_iter!(0..ncol).map(|j| cpp_matrix.entry(i, j)).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        debug_mem(format!("GetMatrixElement row={nrow}, col={ncol}"));
        DCRTPolyMatrix::from_poly_vec(params, matrix_inner)
    }

    pub(crate) fn gadget_vector(params: &DCRTPolyParams) -> DCRTPolyMatrix {
        let base = 1 << params.base_bits();
        let g_vec_cpp = DCRTPolyGadgetVector(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            params.modulus_digits(),
            base,
        );
        DCRTPolyMatrix::from_cpp_matrix_ptr(params, &CppMatrix::new(params, g_vec_cpp))
    }

    fn dcrt_decompose_poly(&self, poly: &DCRTPoly, base_bits: u32) -> Vec<DCRTPoly> {
        let decomposed = poly.get_poly().Decompose(base_bits);
        let cpp_decomposed = CppMatrix::new(&self.params, decomposed);
        parallel_iter!(0..cpp_decomposed.ncol()).map(|idx| cpp_decomposed.entry(0, idx)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{DCRTPolyParams, DCRTPolyUniformSampler, FinRingElem},
        sampler::{DistType, PolyUniformSampler},
    };
    use num_bigint::BigUint;
    use std::sync::Arc;

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
    fn test_matrix_decompose_with_base8() {
        let params = DCRTPolyParams::new(4, 2, 17, 3);
        let digits_length = params.modulus_digits();

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
        assert_eq!(gadget_matrix.size().1, 2 * digits_length);

        let decomposed = matrix.decompose();
        assert_eq!(decomposed.size().0, 2 * digits_length);
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

    #[test]
    fn test_matrix_compact_bytes() {
        let params = DCRTPolyParams::new(8192, 7, 51, 17);
        let sampler = DCRTPolyUniformSampler::new();

        let dists = [DistType::BitDist, DistType::FinRingDist];
        for dist in dists {
            // todo: interesting finding. if its more square shape (e.g (50,50)more than 2m> (100,1)
            // - total 37s) slower
            let ncol = 100;
            let nrow = 1;

            // Create a random matrix
            let matrix = sampler.sample_uniform(&params, nrow, ncol, dist);

            // Convert to compact bytes
            let start_serialize = std::time::Instant::now();
            let compact_bytes = matrix.to_compact_bytes();
            let serialize_time = start_serialize.elapsed();
            println!(
                "to_compact_bytes took: {:?}, bytes_length={}",
                serialize_time,
                compact_bytes.len()
            );

            // Reconstruct from compact bytes
            let start_deserialize = std::time::Instant::now();
            let reconstructed_matrix = DCRTPolyMatrix::from_compact_bytes(&params, &compact_bytes);
            let deserialize_time = start_deserialize.elapsed();
            println!("from_compact_bytes took: {deserialize_time:?}");

            // Verify the matrices are equal
            assert_eq!(matrix, reconstructed_matrix);
        }
    }
}
