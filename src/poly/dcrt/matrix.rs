use num_bigint::BigInt;

use super::{DCRTPoly, DCRTPolyParams, FinRingElem};
use crate::poly::{Poly, PolyMatrix, PolyParams};
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

unsafe impl Send for DCRTPolyMatrix {}
unsafe impl Sync for DCRTPolyMatrix {}

impl Debug for DCRTPolyMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix")
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("inner", &self.inner)
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
}

impl PolyMatrix for DCRTPolyMatrix {
    type P = DCRTPoly;

    fn from_poly_vec(params: &DCRTPolyParams, vec: Vec<Vec<DCRTPoly>>) -> Self {
        let mut c = vec![vec![DCRTPoly::const_zero(params); vec[0].len()]; vec.len()];
        for (i, row) in vec.iter().enumerate() {
            for (j, element) in row.iter().enumerate() {
                c[i][j] = element.clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow: vec.len(), ncol: vec[0].len() }
    }

    fn entry(&self, i: usize, j: usize) -> &Self::P {
        &self.inner[i][j]
    }

    fn get_row(&self, i: usize) -> Vec<Self::P> {
        self.inner[i].clone()
    }

    fn get_column(&self, j: usize) -> Vec<Self::P> {
        let mut column = Vec::with_capacity(self.nrow);
        for i in 0..self.nrow {
            column.push(self.inner[i][j].clone());
        }
        column
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
        let mut c = vec![
            vec![DCRTPoly::const_zero(&self.params); column_end - column_start];
            row_end - row_start
        ];
        for i in row_start..row_end {
            for j in column_start..column_end {
                c[i - row_start][j - column_start] = self.inner[i][j].clone();
            }
        }
        DCRTPolyMatrix {
            inner: c,
            params: self.params.clone(),
            nrow: row_end - row_start,
            ncol: column_end - column_start,
        }
    }

    fn zero(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        let mut c = vec![vec![DCRTPoly::const_zero(params); ncol]; nrow];
        for i in 0..nrow {
            for j in 0..ncol {
                c[i][j] = DCRTPoly::const_zero(params).clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow, ncol }
    }

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        let nrow = size;
        let ncol = size;
        let mut result = vec![vec![DCRTPoly::const_zero(params); ncol]; nrow];
        let scalar = scalar.unwrap_or_else(|| DCRTPoly::const_one(params));
        for i in 0..size {
            result[i][i] = scalar.clone();
        }
        DCRTPolyMatrix { inner: result, params: params.clone(), nrow, ncol }
    }

    fn transpose(&self) -> Self {
        let nrow = self.ncol;
        let ncol = self.nrow;
        let mut result: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); ncol]; nrow];
        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] = self.inner[j][i].clone();
            }
        }
        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[Self]) -> Self {
        #[cfg(debug_assertions)]
        for (idx, other) in others.iter().enumerate() {
            if self.nrow != other.nrow {
                panic!("Concat error: while the shape of the first matrix is ({0}, {1}), that of the {2}-th matirx is ({3},{4})",self.nrow,self.ncol,idx,other.nrow,other.ncol);
            }
        }
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();
        let mut result: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); ncol]; self.nrow];

        // Copy elements from self
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                result[i][j] = self.inner[i][j].clone();
            }
        }

        // Copy elements from others
        let mut offset = self.ncol;
        for other in others {
            for i in 0..self.nrow {
                for j in 0..other.ncol {
                    result[i][offset + j] = other.inner[i][j].clone();
                }
            }
            offset += other.ncol;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow: self.nrow, ncol }
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[Self]) -> Self {
        #[cfg(debug_assertions)]
        for (idx, other) in others.iter().enumerate() {
            if self.ncol != other.ncol {
                panic!("Concat error: while the shape of the first matrix is ({0}, {1}), that of the {2}-th matirx is ({3},{4})",self.nrow,self.ncol,idx,other.nrow,other.ncol);
            }
        }
        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let mut result: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); self.ncol]; nrow];

        // Copy elements from self
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                result[i][j] = self.inner[i][j].clone();
            }
        }

        // Copy elements from others
        let mut offset = self.nrow;
        for other in others {
            for i in 0..other.nrow {
                for j in 0..other.ncol {
                    result[offset + i][j] = other.inner[i][j].clone();
                }
            }
            offset += other.nrow;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol: self.ncol }
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[Self]) -> Self {
        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();
        let mut result: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); ncol]; nrow];

        // Copy elements from self
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                result[i][j] = self.inner[i][j].clone();
            }
        }

        // Copy elements from others
        let mut row_offset = self.nrow;
        let mut col_offset = self.ncol;
        for other in others {
            for i in 0..other.nrow {
                for j in 0..other.ncol {
                    result[row_offset + i][col_offset + j] = other.inner[i][j].clone();
                }
            }
            row_offset += other.nrow;
            col_offset += other.ncol;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    fn tensor(&self, other: &Self) -> Self {
        let nrow = self.nrow * other.nrow;
        let ncol = self.ncol * other.ncol;
        let mut result: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); ncol]; nrow];

        for i1 in 0..self.nrow {
            for j1 in 0..self.ncol {
                for i2 in 0..other.nrow {
                    for j2 in 0..other.ncol {
                        let i = i1 * other.nrow + i2;
                        let j = j1 * other.ncol + j2;
                        result[i][j] = self.inner[i1][j1].clone() * other.inner[i2][j2].clone();
                    }
                }
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

    fn decompose(&self) -> Self {
        let bit_length = self.params.modulus_bits();
        let new_nrow = self.nrow * bit_length;
        let mut new_inner = vec![vec![DCRTPoly::const_zero(&self.params); self.ncol]; new_nrow];

        for i in 0..self.nrow {
            for j in 0..self.ncol {
                let decomposed = self.inner[i][j].decompose(&self.params);
                for bit in 0..bit_length {
                    new_inner[i * bit_length + bit][j] = decomposed[bit].clone();
                }
            }
        }
        Self { nrow: new_nrow, ncol: self.ncol, inner: new_inner, params: self.params.clone() }
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
}

// ====== Arithmetic ======

impl Add for DCRTPolyMatrix {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self + &rhs
    }
}

// Implement addition of a matrix by a matrix reference
impl<'a> Add<&'a DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn add(self, rhs: &'a DCRTPolyMatrix) -> Self::Output {
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
        let mut c: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); self.ncol]; self.nrow];
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                c[i][j] = -self.inner[i][j].clone();
            }
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

// Implement multiplication of a matrix by a matrix reference
impl<'a> Mul<&'a DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: &'a DCRTPolyMatrix) -> Self::Output {
        let nrow = self.nrow;
        let ncol = rhs.ncol;
        #[cfg(debug_assertions)]
        if rhs.nrow != self.ncol {
            panic!(
                "Multiplication condition failed: rhs.nrow ({}) must equal self.ncol ({})",
                rhs.nrow, self.ncol
            );
        }
        let common = self.ncol;
        let mut c: Vec<Vec<DCRTPoly>> = vec![vec![DCRTPoly::const_zero(&self.params); ncol]; nrow];
        for i in 0..nrow {
            for j in 0..ncol {
                for k in 0..common {
                    c[i][j] += self.inner[i][k].clone() * rhs.inner[k][j].clone();
                }
            }
        }
        DCRTPolyMatrix { inner: c, params: self.params, nrow, ncol }
    }
}

// Implement multiplication of a matrix by a polynomial
impl Mul<DCRTPoly> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: DCRTPoly) -> Self::Output {
        self * &rhs
    }
}

// Implement multiplication of a matrix by a polynomial reference
impl Mul<&DCRTPoly> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: &DCRTPoly) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let mut result: Vec<Vec<DCRTPoly>> = self.inner;

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] *= rhs.clone();
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params, nrow, ncol }
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
impl<'a> Sub<&'a DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn sub(self, rhs: &'a DCRTPolyMatrix) -> Self::Output {
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
mod tests {
    use std::sync::Arc;

    use num_bigint::BigUint;

    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;

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
        let decomposed = matrix.decompose();
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
        let prod = matrix1.clone() * &identity;
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
        let col_concat = matrix1.clone().concat_columns(&[matrix2.clone()]);
        assert_eq!(col_concat.row_size(), 2);
        assert_eq!(col_concat.col_size(), 4);

        // Test row concatenation
        let row_concat = matrix1.clone().concat_rows(&[matrix2.clone()]);
        assert_eq!(row_concat.row_size(), 4);
        assert_eq!(row_concat.col_size(), 2);

        // Test diagonal concatenation
        let diag_concat = matrix1.concat_diag(&[matrix2]);
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
        let modulus = params.modulus();

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
    fn test_matrix_addition_mismatch() {
        let params = DCRTPolyParams::default();
        let matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        let matrix2 = DCRTPolyMatrix::zero(&params, 2, 3);
        let _sum = matrix1 + matrix2;
    }

    #[test]
    #[should_panic(expected = "Multiplication condition failed")]
    fn test_matrix_multiplication_mismatch() {
        let params = DCRTPolyParams::default();
        let matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
        let matrix2 = DCRTPolyMatrix::zero(&params, 3, 2);
        let _prod = matrix1 * matrix2;
    }
}
