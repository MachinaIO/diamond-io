use num_bigint::{BigInt, BigUint};
use rayon::iter::{
    IndexedParallelIterator, IntoParallelIterator, IntoParallelRefIterator,
    IntoParallelRefMutIterator, ParallelIterator,
};

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
            .field("params", &self.params)
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
        let nrow = vec.len();
        let ncol = vec[0].len();
        let c: Vec<Vec<DCRTPoly>> = vec
            .into_par_iter()
            .map(|row| row.into_par_iter().map(|element| element).collect())
            .collect();
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow, ncol }
    }

    fn entry(&self, i: usize, j: usize) -> &Self::P {
        &self.inner[i][j]
    }

    fn get_row(&self, i: usize) -> Vec<Self::P> {
        self.inner[i].clone()
    }

    fn get_column(&self, j: usize) -> Vec<Self::P> {
        (0..self.nrow).into_par_iter().map(|i| self.inner[i][j].clone()).collect()
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
        let c: Vec<Vec<DCRTPoly>> = (row_start..row_end)
            .into_par_iter()
            .map(|i| {
                (column_start..column_end)
                    .into_par_iter()
                    .map(|j| self.inner[i][j].clone())
                    .collect()
            })
            .collect();
        DCRTPolyMatrix {
            inner: c,
            params: self.params.clone(),
            nrow: row_end - row_start,
            ncol: column_end - column_start,
        }
    }

    fn zero(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        let c: Vec<Vec<DCRTPoly>> = (0..nrow)
            .into_par_iter()
            .map(|_| (0..ncol).into_par_iter().map(|_| DCRTPoly::const_zero(params)).collect())
            .collect();
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow, ncol }
    }

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        let scalar = scalar.unwrap_or_else(|| DCRTPoly::const_one(params));
        let result: Vec<Vec<DCRTPoly>> = (0..size)
            .into_par_iter()
            .map(|i| {
                let mut row = vec![DCRTPoly::const_zero(params); size];
                row[i] = scalar.clone();
                row
            })
            .collect();
        DCRTPolyMatrix { inner: result, params: params.clone(), nrow: size, ncol: size }
    }

    fn transpose(&self) -> Self {
        let nrow = self.ncol;
        let ncol = self.nrow;
        let result: Vec<Vec<DCRTPoly>> = (0..nrow)
            .into_par_iter()
            .map(|i| (0..ncol).into_par_iter().map(|j| self.inner[j][i].clone()).collect())
            .collect();
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

        let result: Vec<Vec<DCRTPoly>> = (0..nrow)
            .into_par_iter()
            .map(|i| {
                (0..ncol)
                    .into_par_iter()
                    .map(|j| {
                        let i1 = i / other.nrow;
                        let i2 = i % other.nrow;
                        let j1 = j / other.ncol;
                        let j2 = j % other.ncol;
                        self.inner[i1][j1].clone() * other.inner[i2][j2].clone()
                    })
                    .collect()
            })
            .collect();

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self {
        let bit_length = params.modulus_bits();
        let poly_vec: Vec<DCRTPoly> = (0u32..bit_length as u32)
            .into_par_iter()
            .map(|i| {
                let value = BigInt::from(2).pow(i);
                DCRTPoly::from_const(params, &FinRingElem::new(value, params.modulus()))
            })
            .collect();
        let gadget_vector = Self::from_poly_vec(params, vec![poly_vec]);
        let identity = DCRTPolyMatrix::identity(params, size, None);
        identity.tensor(&gadget_vector)
    }

    fn decompose(&self) -> Self {
        let bit_length = self.params.modulus_bits();
        let new_nrow = self.nrow * bit_length;

        let new_inner: Vec<Vec<DCRTPoly>> = (0..new_nrow)
            .into_par_iter()
            .map(|new_i| {
                let i = new_i / bit_length;
                let bit = new_i % bit_length;
                (0..self.ncol)
                    .into_par_iter()
                    .map(|j| {
                        let coeffs = self.inner[i][j].coeffs();
                        let bit_coeffs: Vec<_> = coeffs
                            .par_iter()
                            .map(|coeff_val| {
                                let val = (coeff_val.value() >> bit) & BigUint::from(1u32);
                                FinRingElem::new(val, self.params.modulus())
                            })
                            .collect();
                        DCRTPoly::from_coeffs(&self.params, &bit_coeffs)
                    })
                    .collect()
            })
            .collect();

        Self { nrow: new_nrow, ncol: self.ncol, inner: new_inner, params: self.params.clone() }
    }
}

// ====== Arithmetic ======

impl Add for DCRTPolyMatrix {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self + &rhs
    }
}

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

        let result: Vec<Vec<DCRTPoly>> = self
            .inner
            .into_par_iter()
            .zip(rhs.inner.par_iter())
            .map(|(self_row, rhs_row)| {
                self_row
                    .into_par_iter()
                    .zip(rhs_row.par_iter())
                    .map(|(a, b)| a + b.clone())
                    .collect()
            })
            .collect();

        Self { inner: result, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}

impl Neg for DCRTPolyMatrix {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let c: Vec<Vec<DCRTPoly>> = self
            .inner
            .into_par_iter()
            .map(|row| row.into_par_iter().map(|element| -element).collect())
            .collect();

        DCRTPolyMatrix { inner: c, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}

impl Mul for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
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

        let c: Vec<Vec<DCRTPoly>> = (0..nrow)
            .into_par_iter()
            .map(|i| {
                (0..ncol)
                    .into_par_iter()
                    .map(|j| {
                        (0..common)
                            .map(|k| self.inner[i][k].clone() * rhs.inner[k][j].clone())
                            .reduce(|acc, x| acc + x)
                            .expect("Cannot accumulate elements while multiplication")
                    })
                    .collect()
            })
            .collect();

        Self { inner: c, params: self.params, nrow, ncol }
    }
}

impl Mul<DCRTPoly> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(self, rhs: DCRTPoly) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&DCRTPoly> for DCRTPolyMatrix {
    type Output = Self;

    fn mul(mut self, rhs: &DCRTPoly) -> Self::Output {
        self.inner.par_iter_mut().for_each(|row| {
            row.par_iter_mut().for_each(|element| {
                *element *= rhs.clone();
            });
        });

        Self { inner: self.inner, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}
// Implement subtraction for matrices
impl Sub for DCRTPolyMatrix {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<'a> Sub<&'a DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = Self;

    fn sub(mut self, rhs: &'a DCRTPolyMatrix) -> Self::Output {
        #[cfg(debug_assertions)]
        if self.nrow != rhs.nrow || self.ncol != rhs.ncol {
            panic!(
                "Subtraction requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
                self.nrow, self.ncol, rhs.nrow, rhs.ncol
            );
        }

        self.inner.par_iter_mut().enumerate().for_each(|(i, row)| {
            row.iter_mut().enumerate().for_each(|(j, element)| {
                *element -= rhs.inner[i][j].clone();
            });
        });

        DCRTPolyMatrix { inner: self.inner, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}

#[cfg(test)]
mod tests {
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
