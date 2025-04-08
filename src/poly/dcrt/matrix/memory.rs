use crate::{
    parallel_iter,
    poly::{
        dcrt::{cpp_matrix::CppMatrix, DCRTPoly, DCRTPolyParams, FinRingElem},
        Poly, PolyMatrix, PolyParams,
    },
};
use num_bigint::BigInt;
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
};

#[derive(Clone)]
pub struct DCRTPolyMatrix {
    inner: Vec<Vec<DCRTPoly>>,
    params: DCRTPolyParams,
    nrow: usize,
    ncol: usize,
}

impl Debug for DCRTPolyMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix")
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("params", &self.params)
            .field("inner", &self.inner)
            .finish()
    }
}

impl PartialEq for DCRTPolyMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner && self.nrow == other.nrow && self.ncol == other.ncol
    }
}

impl Eq for DCRTPolyMatrix {}

// Add getter methods for inner and params
impl DCRTPolyMatrix {
    pub fn inner(&self) -> &Vec<Vec<DCRTPoly>> {
        &self.inner
    }
    pub fn params(&self) -> &DCRTPolyParams {
        &self.params
    }

    pub(crate) fn dcrt_decompose_poly(poly: &DCRTPoly, base_bits: u32) -> Vec<DCRTPoly> {
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

    fn decompose(&self, base_bits: Option<u32>) -> Self {
        let bit_length = self.params.modulus_bits();
        let base_bits = base_bits.unwrap_or(self.params.base_bits());

        Self {
            nrow: self.nrow * bit_length,
            ncol: self.ncol,
            inner: parallel_iter!(0..self.nrow)
                .flat_map(|i| {
                    let decompositions: Vec<_> = parallel_iter!(0..self.ncol)
                        .map(|j| Self::dcrt_decompose_poly(&self.inner[i][j], base_bits))
                        .collect();
                    let mut decomposed_rows = vec![Vec::with_capacity(self.ncol); bit_length];
                    for col_decomp in decompositions {
                        for bit in 0..bit_length {
                            decomposed_rows[bit].push(col_decomp[bit].clone());
                        }
                    }
                    decomposed_rows
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
        assert_eq!(self.ncol, other.nrow * identity_size * self.params.modulus_bits());
        let slice_width = other.ncol;
        let mut output = vec![DCRTPolyMatrix::zero(&self.params, 0, 0); other.ncol * identity_size];

        for j in 0..other.ncol {
            let jth_col_m_decompose = other.get_column_matrix_decompose(j, None);
            for i in 0..identity_size {
                let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
                output[i * other.ncol + j] = slice * &jth_col_m_decompose;
            }
        }

        output[0].clone().concat_columns(&output[1..].iter().collect::<Vec<_>>())
    }

    fn get_column_matrix_decompose(&self, j: usize, base_bits: Option<u32>) -> Self {
        DCRTPolyMatrix::from_poly_vec(
            &self.params,
            self.get_column(j).into_iter().map(|poly| vec![poly]).collect(),
        )
        .decompose(base_bits)
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

// Implement subtraction of a matrix by a matrix reference
impl Sub<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

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
        let mut result = self.inner;

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] -= rhs.inner[i][j].clone();
            }
        }

        Self { inner: result, params: self.params, ncol, nrow }
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
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
    fn test_decompose() {
        let params = DCRTPolyParams::default();
        let bit_length = params.modulus_bits();

        // Create a simple 2x8 matrix with some non-zero values
        let mut matrix = DCRTPolyMatrix::zero(&params, 2, 8);
        assert_eq!(matrix.row_size(), 2);
        assert_eq!(matrix.col_size(), 8);
        let value = FinRingElem::new(5u32, params.modulus());
        matrix.inner[0][0] = DCRTPoly::from_const(&params, &value);
        matrix.inner[1][1] = DCRTPoly::from_const(&params, &value);
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, 2);
        assert_eq!(gadget_matrix.row_size(), 2);
        assert_eq!(gadget_matrix.col_size(), 2 * bit_length);
        let decomposed = matrix.decompose(None);
        assert_eq!(decomposed.row_size(), 2 * bit_length);
        assert_eq!(decomposed.col_size(), 8);

        let expected_matrix = gadget_matrix * decomposed;
        assert_eq!(expected_matrix.row_size(), 2);
        assert_eq!(expected_matrix.col_size(), 8);
        assert_eq!(matrix, expected_matrix);
    }

    #[test]
    fn test_matrix_basic_operations() {
        let params = DCRTPolyParams::default();

        // Test zero and identity matrices
        let zero = DCRTPolyMatrix::zero(&params, 2, 2);
        let identity = DCRTPolyMatrix::identity(&params, 2, None);

        // Test matrix creation and equality
        let mut matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        let value = FinRingElem::new(5u32, params.modulus());
        matrix1.inner[0][0] = DCRTPoly::from_const(&params, &value);
        matrix1.inner[1][1] = DCRTPoly::from_const(&params, &value);

        let matrix2 = matrix1.clone();
        assert_eq!(matrix1, matrix2);

        // Test addition
        let sum = matrix1.clone() + &matrix2;
        let value_10 = FinRingElem::new(10u32, params.modulus());
        let _expected_sum = DCRTPolyMatrix::zero(&params, 2, 2);
        assert_eq!(sum.entry(0, 0).coeffs()[0], value_10);

        // Test subtraction
        let diff = matrix1.clone() - &matrix2;
        assert_eq!(diff, zero);

        // Test multiplication
        let prod = &matrix1 * &identity;
        assert_eq!(prod, matrix1);
    }

    #[test]
    fn test_matrix_concatenation() {
        let params = DCRTPolyParams::default();
        let value = FinRingElem::new(5u32, params.modulus());

        let mut matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        matrix1.inner[0][0] = DCRTPoly::from_const(&params, &value);

        let mut matrix2 = DCRTPolyMatrix::zero(&params, 2, 2);
        matrix2.inner[1][1] = DCRTPoly::from_const(&params, &value);

        // Test column concatenation
        let col_concat = matrix1.concat_columns(&[&matrix2]);
        assert_eq!(col_concat.row_size(), 2);
        assert_eq!(col_concat.col_size(), 4);

        // Test row concatenation
        let row_concat = matrix1.concat_rows(&[&matrix2]);
        assert_eq!(row_concat.row_size(), 4);
        assert_eq!(row_concat.col_size(), 2);

        // Test diagonal concatenation
        let diag_concat = matrix1.concat_diag(&[&matrix2]);
        assert_eq!(diag_concat.row_size(), 4);
        assert_eq!(diag_concat.col_size(), 4);
    }

    #[test]
    fn test_matrix_tensor_product() {
        let params = DCRTPolyParams::default();
        let value = FinRingElem::new(5u32, params.modulus());

        let mut matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        matrix1.inner[0][0] = DCRTPoly::from_const(&params, &value);

        let mut matrix2 = DCRTPolyMatrix::zero(&params, 2, 2);
        matrix2.inner[0][0] = DCRTPoly::from_const(&params, &value);

        let tensor = matrix1.tensor(&matrix2);
        assert_eq!(tensor.row_size(), 4);
        assert_eq!(tensor.col_size(), 4);

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

    // #[test]
    // fn test_mul_tensor_identity() {
    //     let params = DCRTPolyParams::default();
    //     let sampler = DCRTPolyUniformSampler::new();

    //     // Create matrix S (2x12)
    //     let s = sampler.sample_uniform(&params, 2, 12,
    // crate::poly::sampler::DistType::FinRingDist);

    //     // Create 'other' matrix (3x3)
    //     let other =
    //         sampler.sample_uniform(&params, 3, 3, crate::poly::sampler::DistType::FinRingDist);

    //     // Perform S * (I_4 ⊗ other)
    //     let result = s.mul_tensor_identity(&other, 4);

    //     // Check dimensions
    //     assert_eq!(result.row_size(), 2);
    //     assert_eq!(result.col_size(), 12);

    //     let identity = DCRTPolyMatrix::identity(&params, 4, None);

    //     // Check result
    //     let expected_result = s * (identity.tensor(&other));

    //     assert_eq!(expected_result.row_size(), 2);
    //     assert_eq!(expected_result.col_size(), 12);
    //     assert_eq!(result, expected_result)
    // }

    // #[test]
    // fn test_mul_tensor_identity_decompose_naive() {
    //     let params = DCRTPolyParams::default();
    //     let sampler = DCRTPolyUniformSampler::new();

    //     // Create matrix S (2x2516)
    //     let s =
    //         sampler.sample_uniform(&params, 2, 2516,
    // crate::poly::sampler::DistType::FinRingDist);

    //     // Create 'other' matrix (2x68)
    //     let other =
    //         sampler.sample_uniform(&params, 2, 68, crate::poly::sampler::DistType::FinRingDist);

    //     // Decompose 'other' matrix
    //     let other_decompose = other.decompose();

    //     // Perform S * (I_37 ⊗ G^-1(other))
    //     let result: DCRTPolyMatrix = s.mul_tensor_identity(&other_decompose, 37);

    //     // Check dimensions
    //     assert_eq!(result.row_size(), 2);
    //     assert_eq!(result.col_size(), 2516);

    //     let identity = DCRTPolyMatrix::identity(&params, 37, None);

    //     // Check result
    //     let expected_result = s * (identity.tensor(&other_decompose));

    //     assert_eq!(expected_result.row_size(), 2);
    //     assert_eq!(expected_result.col_size(), 2516);
    //     assert_eq!(result, expected_result)
    // }

    // #[test]
    // fn test_mul_tensor_identity_decompose_optimal() {
    //     let params = DCRTPolyParams::default();
    //     let sampler = DCRTPolyUniformSampler::new();

    //     // Create matrix S (2x2516)
    //     let s =
    //         sampler.sample_uniform(&params, 2, 2516,
    // crate::poly::sampler::DistType::FinRingDist);

    //     // Create 'other' matrix (2x68)
    //     let other =
    //         sampler.sample_uniform(&params, 2, 68, crate::poly::sampler::DistType::FinRingDist);

    //     // Perform S * (I_37 ⊗ G^-1(other))
    //     let result: DCRTPolyMatrix = s.mul_tensor_identity_decompose(&other, 37);

    //     // // Check dimensions
    //     assert_eq!(result.row_size(), 2);
    //     assert_eq!(result.col_size(), 2516);

    //     let identity = DCRTPolyMatrix::identity(&params, 37, None);

    //     // Check result
    //     let expected_result = s.clone() * (identity.tensor(&other.decompose()));
    //     let expected_result_2 = s.mul_tensor_identity(&other.decompose(), 37);

    //     assert_eq!(expected_result.row_size(), 2);
    //     assert_eq!(expected_result.col_size(), 2516);

    //     assert_eq!(expected_result_2.row_size(), 2);
    //     assert_eq!(expected_result_2.col_size(), 2516);

    //     assert_eq!(result, expected_result);
    //     assert_eq!(result, expected_result_2);
    // }
}
