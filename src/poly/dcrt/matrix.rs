use num_bigint::{BigInt, BigUint};

use super::{DCRTPoly, DCRTPolyParams, FinRing};
use crate::{
    poly::{Poly, PolyMatrix, PolyParams},
    utils::ceil_log2,
};
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

    /// Gadget vector g = (2^0, 2^1, ..., 2^{log(q)-1})
    /// where g ∈ Z_q^{log(q)}
    fn gadget_vector(params: &<Self::P as Poly>::Params) -> Self {
        let q = params.modulus();
        let size = ceil_log2(&q);
        let mut poly_vec = Vec::with_capacity(size);
        for i in 0..size {
            let value = BigInt::from(2).pow(i.try_into().unwrap());
            let fe: FinRing = FinRing::new(value, q.clone().into());
            poly_vec.push(DCRTPoly::from_const(params, &fe));
        }
        Self::from_poly_vec(params, vec![poly_vec])
    }

    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self {
        let identity = DCRTPolyMatrix::identity(params, size, None);
        let gadget_vector = Self::gadget_vector(params);
        identity.tensor(&gadget_vector.transpose())
    }

    fn decompose(&self) -> Self {
        let q = self.params.modulus();
        let bit_length = ceil_log2(&q);
        let new_ncol = self.ncol * bit_length;
        let mut new_inner = Vec::with_capacity(self.nrow);

        for i in 0..self.nrow {
            let mut new_row = Vec::with_capacity(new_ncol);
            for _ in 0..new_ncol {
                new_row.push(DCRTPoly::const_zero(&self.params));
            }
            for j in 0..self.ncol {
                let c_ij = &self.inner[i][j];
                let coeffs = c_ij.coeffs();
                let coeff_len = coeffs.len();
                for bit in 0..bit_length {
                    let mut bit_coeffs = Vec::with_capacity(coeff_len);
                    for coeff_val in coeffs.clone() {
                        // bit_value in {0, 1}
                        let val = (coeff_val.value() >> bit) & BigUint::from(1u32);
                        let elem = FinRing::new(val.clone(), self.params.modulus().into());
                        bit_coeffs.push(elem);
                    }
                    let bit_poly = DCRTPoly::from_coeffs(&self.params, &bit_coeffs);
                    new_row[j * bit_length + bit] = bit_poly;
                }
            }
            new_inner.push(new_row);
        }
        Self { nrow: self.nrow, ncol: new_ncol, inner: new_inner, params: self.params.clone() }
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
    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;

    #[test]
    fn test_gadget_vector() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let gadget_vector = DCRTPolyMatrix::gadget_vector(&params);
        assert_eq!(gadget_vector.row_size(), 1);
        assert_eq!(gadget_vector.col_size(), ceil_log2(&params.modulus()) + 1);
    }

    #[test]
    fn test_gadget_matrix() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let size = 3;
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, size);
        assert_eq!(gadget_matrix.row_size(), size * (ceil_log2(&params.modulus()) + 1));
        assert_eq!(gadget_matrix.col_size(), size);
    }

    #[test]
    fn test_decompose() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let bit_length = ceil_log2(&params.modulus());

        // Create a simple 2x2 matrix with some non-zero values
        let mut matrix = DCRTPolyMatrix::zero(&params, 2, 2);
        let value = FinRing::new(5u32, params.modulus().into());
        matrix.inner[0][0] = DCRTPoly::from_const(&params, &value);
        matrix.inner[1][1] = DCRTPoly::from_const(&params, &value);

        let decomposed = matrix.decompose();

        // Check dimensions
        assert_eq!(decomposed.row_size(), 2);
        assert_eq!(decomposed.col_size(), 2 * bit_length);
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, 2);
        let expected_matrix = decomposed * gadget_matrix;
        assert_eq!(matrix, expected_matrix);
    }
}
