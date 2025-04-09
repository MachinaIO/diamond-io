use crate::{
    parallel_iter,
    poly::{
        dcrt::{cpp_matrix::CppMatrix, DCRTPoly, DCRTPolyParams},
        MatrixElem, Poly, PolyMatrix, PolyParams,
    },
    utils::debug_mem,
};
use itertools::Itertools;
use openfhe::ffi::{DCRTPolyGadgetVector, MatrixGen, SetMatrixElement};
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Range, Sub},
};

#[derive(Clone)]
pub struct DCRTPolyMatrix<T: MatrixElem> {
    inner: Vec<Vec<T>>,
    pub params: DCRTPolyParams,
    nrow: usize,
    ncol: usize,
}

impl<T: MatrixElem> Debug for DCRTPolyMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix")
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("params", &self.params)
            .field("inner", &self.inner)
            .finish()
    }
}

impl<T: MatrixElem> PartialEq for DCRTPolyMatrix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner && self.nrow == other.nrow && self.ncol == other.ncol
    }
}

impl<T: MatrixElem> Eq for DCRTPolyMatrix<T> {}

// Add getter methods for inner and params
impl<T: MatrixElem> DCRTPolyMatrix<T> {
    pub fn inner(&self) -> &Vec<Vec<T>> {
        &self.inner
    }

    pub fn new_empty(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        let inner = vec![vec![T::new_empty(); ncol]; nrow];
        Self { inner, params: params.clone(), nrow, ncol }
    }

    pub fn replace_entries<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<DCRTPoly>> + Send + Sync,
    {
        let new_entries = f(rows.clone(), cols.clone());
        assert_eq!(
            new_entries.len(),
            rows.end - rows.start,
            "Returned number of rows does not match the provided row range"
        );
        for (i, row_vec) in new_entries.iter().enumerate() {
            assert_eq!(
                row_vec.len(),
                cols.end - cols.start,
                "Returned number of columns does not match the provided column range at row {}",
                i
            );
        }
        for (i, row_entries) in new_entries.into_iter().enumerate() {
            let row_idx = rows.start + i;
            for (j, entry) in row_entries.into_iter().enumerate() {
                let col_idx = cols.start + j;
                self.inner[row_idx][col_idx] = entry;
            }
        }
    }

    pub fn block_entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<DCRTPoly>> {
        rows.into_par_iter()
            .map(|i| {
                // Use slice indexing on the row. The slice is converted into a Vec by `.to_vec()`.
                self.inner[i][cols.clone()].to_vec()
            })
            .collect()
    }

    pub(crate) fn to_cpp_matrix_ptr(&self) -> CppMatrix {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let params = &self.params;
        let mut matrix_ptr =
            MatrixGen(params.ring_dimension(), params.crt_depth(), params.crt_bits(), nrow, ncol);
        debug_mem(format!("matrix_ptr MatrixGen row={}, col={}", nrow, ncol));
        for i in 0..nrow {
            for j in 0..ncol {
                SetMatrixElement(matrix_ptr.as_mut().unwrap(), i, j, self.entry(i, j).get_poly());
            }
        }
        debug_mem(format!("SetMatrixElement row={}, col={}", nrow, ncol));
        CppMatrix::new(matrix_ptr)
    }

    pub(crate) fn from_cpp_matrix_ptr(params: &DCRTPolyParams, cpp_matrix: &CppMatrix) -> Self {
        let nrow = cpp_matrix.nrow();
        let ncol = cpp_matrix.ncol();
        let matrix_inner = parallel_iter!(0..nrow)
            .map(|i| parallel_iter!(0..ncol).map(|j| cpp_matrix.entry(i, j)).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        debug_mem(format!("GetMatrixElement row={}, col={}", nrow, ncol));
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
        DCRTPolyMatrix::from_cpp_matrix_ptr(params, &CppMatrix::new(g_vec_cpp))
    }

    fn dcrt_decompose_poly(poly: &DCRTPoly, base_bits: u32) -> Vec<DCRTPoly> {
        let decomposed = poly.get_poly().Decompose(base_bits);
        let cpp_decomposed = CppMatrix::new(decomposed);
        parallel_iter!(0..cpp_decomposed.ncol()).map(|idx| cpp_decomposed.entry(0, idx)).collect()
    }
}

impl PolyMatrix for DCRTPolyMatrix {
    type P = DCRTPoly;

    fn from_poly_vec(params: &DCRTPolyParams, vec: Vec<Vec<DCRTPoly>>) -> Self {
        let nrow = vec.len();
        let ncol = vec[0].len();
        DCRTPolyMatrix { inner: vec, params: params.clone(), nrow, ncol }
    }

    fn entry(&self, i: usize, j: usize) -> Self::P {
        self.inner[i][j].clone()
    }

    fn get_row(&self, i: usize) -> Vec<Self::P> {
        self.inner[i].clone()
    }

    fn get_column(&self, j: usize) -> Vec<Self::P> {
        self.inner.iter().map(|row| row[j].clone()).collect()
    }

    fn size(&self) -> (usize, usize) {
        (self.nrow, self.ncol)
    }

    fn row_size(&self) -> usize {
        self.nrow
    }

    fn col_size(&self) -> usize {
        self.ncol
    }

    fn slice(
        &self,
        row_start: usize,
        row_end: usize,
        column_start: usize,
        column_end: usize,
    ) -> Self {
        let nrow = row_end - row_start;
        let ncol = column_end - column_start;

        let mut c = Vec::with_capacity(nrow);
        for i in row_start..row_end {
            let mut row = Vec::with_capacity(ncol);
            for j in column_start..column_end {
                row.push(self.inner[i][j].clone());
            }
            c.push(row);
        }

        DCRTPolyMatrix { inner: c, params: self.params.clone(), nrow, ncol }
    }

    fn zero(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        let mut c = Vec::with_capacity(nrow);
        let zero_elem = DCRTPoly::const_zero(params);
        for _ in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for _ in 0..ncol {
                row.push(zero_elem.clone());
            }
            c.push(row);
        }
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow, ncol }
    }

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        let nrow = size;
        let ncol = size;
        let scalar = scalar.unwrap_or_else(|| DCRTPoly::const_one(params));
        let zero_elem = DCRTPoly::const_zero(params);
        let mut result = Vec::with_capacity(nrow);

        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                if i == j {
                    row.push(scalar.clone());
                } else {
                    row.push(zero_elem.clone());
                }
            }
            result.push(row);
        }

        DCRTPolyMatrix { inner: result, params: params.clone(), nrow, ncol }
    }

    fn transpose(&self) -> Self {
        let nrow = self.ncol;
        let ncol = self.nrow;
        let mut result = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                row.push(self.inner[j][i].clone());
            }
            result.push(row);
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[&Self]) -> Self {
        #[cfg(debug_assertions)]
        for (idx, other) in others.iter().enumerate() {
            if self.nrow != other.nrow {
                panic!(
                    "Concat error: while the shape of the first matrix is ({}, {}), that of the {}-th matrix is ({},{})",
                    self.nrow, self.ncol, idx, other.nrow, other.ncol
                );
            }
        }
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();
        let result: Vec<Vec<DCRTPoly>> = parallel_iter!(0..self.nrow)
            .map(|i| {
                let mut row: Vec<&DCRTPoly> = Vec::with_capacity(ncol);
                row.extend(self.inner[i].iter());
                for other in others {
                    row.extend(other.inner[i].iter());
                }
                row.into_iter().cloned().collect()
            })
            .collect();

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow: self.nrow, ncol }
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[&Self]) -> Self {
        #[cfg(debug_assertions)]
        for (idx, other) in others.iter().enumerate() {
            if self.ncol != other.ncol {
                panic!(
                    "Concat error: while the shape of the first matrix is ({}, {}), that of the {}-th matrix is ({}, {})",
                    self.nrow, self.ncol, idx, other.nrow, other.ncol
                );
            }
        }

        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let mut result = Vec::with_capacity(nrow);
        for i in 0..self.nrow {
            let mut row = Vec::with_capacity(self.ncol);
            for j in 0..self.ncol {
                row.push(self.inner[i][j].clone());
            }
            result.push(row);
        }
        for other in others {
            for i in 0..other.nrow {
                let mut row = Vec::with_capacity(self.ncol);
                for j in 0..other.ncol {
                    row.push(other.inner[i][j].clone());
                }
                result.push(row);
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol: self.ncol }
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[&Self]) -> Self {
        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();

        let zero_elem = DCRTPoly::const_zero(&self.params);

        let mut result: Vec<Vec<DCRTPoly>> = Vec::with_capacity(nrow);

        // First part of the matrix (self)
        for i in 0..self.nrow {
            let mut row = Vec::with_capacity(ncol);
            row.extend(self.inner[i].iter().cloned());
            row.extend(std::iter::repeat_n(zero_elem.clone(), ncol - self.ncol));
            result.push(row);
        }

        let mut col_offset = self.ncol;
        // Concatenating with others in parallel
        for other in others.iter() {
            result.extend(
                parallel_iter!(0..other.nrow)
                    .map(|i| {
                        let mut row = Vec::with_capacity(ncol);
                        row.extend(std::iter::repeat_n(zero_elem.clone(), col_offset));
                        row.extend(other.inner[i].iter().cloned());
                        row.extend(std::iter::repeat_n(
                            zero_elem.clone(),
                            ncol - col_offset - other.ncol,
                        ));
                        row
                    })
                    .collect::<Vec<Vec<DCRTPoly>>>()
                    .into_iter(),
            );
            col_offset += other.ncol;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    fn tensor(&self, other: &Self) -> Self {
        let nrow = self.nrow * other.nrow;
        let ncol = self.ncol * other.ncol;
        let mut result: Vec<Vec<DCRTPoly>> = Vec::with_capacity(nrow);
        for i1 in 0..self.nrow {
            for i2 in 0..other.nrow {
                let mut row = Vec::with_capacity(ncol);
                for j1 in 0..self.ncol {
                    for j2 in 0..other.ncol {
                        row.push(&self.inner[i1][j1] * &other.inner[i2][j2]);
                    }
                }
                result.push(row);
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
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
        let decomposed_entries: Vec<Vec<Vec<DCRTPoly>>> = parallel_iter!(0..self.nrow)
            .map(|i| {
                (0..self.ncol)
                    .map(|j| Self::dcrt_decompose_poly(&self.entry(i, j), base_bits))
                    .collect()
            })
            .collect();
        Self {
            nrow: self.nrow * log_base_q,
            ncol: self.ncol,
            inner: parallel_iter!(0..self.nrow * log_base_q)
                .map(|idx| {
                    let i = idx / log_base_q;
                    let k = idx % log_base_q;

                    (0..self.ncol).map(|j| decomposed_entries[i][j][k].clone()).collect()
                })
                .collect(),
            params: self.params.clone(),
        }
    }

    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self {
        let mut new_inner = self.clone().inner;
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                new_inner[i][j] =
                    self.inner[i][j].modulus_switch(&self.params, new_modulus.clone());
            }
        }
        Self { inner: new_inner, params: self.params.clone(), nrow: self.nrow, ncol: self.ncol }
    }

    fn mul_tensor_identity(&self, other: &Self, identity_size: usize) -> Self {
        assert_eq!(self.ncol, other.nrow * identity_size);
        let slice_width = other.nrow;
        let mut slice_results = Vec::with_capacity(identity_size);
        for i in 0..identity_size {
            let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
            slice_results.push(slice * other);
        }
        slice_results[0].clone().concat_columns(&slice_results[1..].iter().collect::<Vec<_>>())
    }

    fn mul_tensor_identity_decompose(&self, other: &Self, identity_size: usize) -> Self {
        let log_base_q = self.params.modulus_digits();
        debug_assert_eq!(self.ncol, other.nrow * identity_size * log_base_q);
        let slice_width = other.nrow * log_base_q;

        let output = (0..identity_size)
            .flat_map(|i| {
                let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
                (0..other.ncol).map(move |j| &slice * &other.get_column_matrix_decompose(j))
            })
            .collect_vec();

        output[0].concat_columns(&output[1..].iter().collect::<Vec<_>>())
    }

    fn get_column_matrix_decompose(&self, j: usize) -> Self {
        DCRTPolyMatrix::from_poly_vec(
            &self.params,
            self.get_column(j).into_iter().map(|poly| vec![poly]).collect(),
        )
        .decompose()
    }
}

// ====== Arithmetic ======

impl Add for DCRTPolyMatrix {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self + &rhs
    }
}

// Implement addition of a matrix by a matrix reference
impl Add<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn add(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        #[cfg(debug_assertions)]
        if self.nrow != rhs.nrow || self.ncol != rhs.ncol {
            panic!(
                "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
                self.nrow, self.ncol, rhs.nrow, rhs.ncol
            );
        }

        let nrow = self.row_size();
        let ncol = self.col_size();
        let mut result = self.inner;

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] += rhs.inner[i][j].clone();
            }
        }

        Self { inner: result, params: self.params, ncol, nrow }
    }
}

impl Neg for DCRTPolyMatrix {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut c = Vec::with_capacity(self.nrow);
        for i in 0..self.nrow {
            let mut row = Vec::with_capacity(self.ncol);
            for j in 0..self.ncol {
                row.push(-self.inner[i][j].clone());
            }
            c.push(row);
        }

        DCRTPolyMatrix { inner: c, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}

impl Mul for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl Mul<&DCRTPolyMatrix> for &DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn mul(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        let nrow = self.nrow;
        let ncol = rhs.ncol;
        #[cfg(debug_assertions)]
        if rhs.nrow != self.ncol {
            panic!(
                "Multiplication condition failed: rhs.nrow ({}) must equal self.ncol ({})",
                rhs.nrow, self.ncol
            );
        }

        DCRTPolyMatrix {
            inner: parallel_iter!(0..nrow)
                .map(|i| {
                    parallel_iter!(0..ncol)
                        .map(|j| {
                            let mut sum = &self.inner[i][0] * &rhs.inner[0][j];
                            for k in 1..self.ncol {
                                sum += &self.inner[i][k] * &rhs.inner[k][j];
                            }
                            sum
                        })
                        .collect()
                })
                .collect(),
            params: self.params.clone(),
            nrow,
            ncol,
        }
    }
}

// Implement multiplication of a matrix by a polynomial
impl Mul<DCRTPoly> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: DCRTPoly) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&DCRTPoly> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: &DCRTPoly) -> Self::Output {
        &self * rhs
    }
}

impl Mul<&DCRTPoly> for &DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn mul(self, rhs: &DCRTPoly) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;

        let result: Vec<Vec<DCRTPoly>> = parallel_iter!(0..nrow)
            .map(|i| parallel_iter!(0..ncol).map(|j| &self.inner[i][j] * rhs).collect())
            .collect();

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }
}

// Implement subtraction for matrices
impl Sub for DCRTPolyMatrix {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl Sub<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn sub(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        &self - rhs
    }
}

// Implement subtraction of a matrix by a matrix reference
impl Sub<&DCRTPolyMatrix> for &DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn sub(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        #[cfg(debug_assertions)]
        if self.nrow != rhs.nrow || self.ncol != rhs.ncol {
            panic!(
                "Subtraction requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
                self.nrow, self.ncol, rhs.nrow, rhs.ncol
            );
        }

        let nrow = self.row_size();
        let ncol = self.col_size();
        let mut result = self.inner.clone();

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] -= rhs.inner[i][j].clone();
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), ncol, nrow }
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use crate::poly::{
        dcrt::{DCRTPolyUniformSampler, FinRingElem},
        sampler::PolyUniformSampler,
    };

    use super::*;
    use num_bigint::BigUint;
    use std::sync::Arc;

    #[test]
    fn test_gadget_matrix() {
        let params = DCRTPolyParams::default();
        let size = 3;
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, size);
        assert_eq!(gadget_matrix.row_size(), size);
        assert_eq!(gadget_matrix.col_size(), size * params.modulus_bits());
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
    fn test_modulus_switch() {
        let params = DCRTPolyParams::default();

        let value00 = FinRingElem::new(1023782870921908217643761278891282178u128, params.modulus());
        let value01 = FinRingElem::new(8179012198875468938912873783289218738u128, params.modulus());
        let value10 = FinRingElem::new(2034903202902173762872163465127672178u128, params.modulus());
        let value11 = FinRingElem::new(1990091289902891278121564387120912660u128, params.modulus());
        let matrix = DCRTPolyMatrix::from_poly_vec(
            &params,
            vec![
                vec![
                    DCRTPoly::from_const(&params, &value00),
                    DCRTPoly::from_const(&params, &value01),
                ],
                vec![
                    DCRTPoly::from_const(&params, &value10),
                    DCRTPoly::from_const(&params, &value11),
                ],
            ],
        );
        let new_modulus = Arc::new(BigUint::from(2u32));
        let switched = matrix.modulus_switch(&new_modulus);
        // although the value becomes less than the new modulus, the set modulus is still the same
        assert_eq!(switched.params().modulus(), params.modulus());
        let new_value00 = value00.modulus_switch(new_modulus.clone());
        let new_value01 = value01.modulus_switch(new_modulus.clone());
        let new_value10 = value10.modulus_switch(new_modulus.clone());
        let new_value11 = value11.modulus_switch(new_modulus.clone());
        let expected = DCRTPolyMatrix::from_poly_vec(
            &params,
            vec![
                vec![
                    DCRTPoly::from_const(&params, &new_value00),
                    DCRTPoly::from_const(&params, &new_value01),
                ],
                vec![
                    DCRTPoly::from_const(&params, &new_value10),
                    DCRTPoly::from_const(&params, &new_value11),
                ],
            ],
        );
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
