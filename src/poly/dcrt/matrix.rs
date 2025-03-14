use super::{DCRTPoly, DCRTPolyParams, FinRingElem};
use crate::{
    parallel_iter,
    poly::{Poly, PolyMatrix, PolyParams},
};
use bytes::Bytes;
use num_bigint::BigInt;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
};
use tracing::info;

#[derive(Clone)]
pub struct DCRTPolyMatrix {
    entries: HashMap<(usize, usize), DCRTPoly>,
    params: DCRTPolyParams,
    nrow: usize,
    ncol: usize,
    zero_elem: DCRTPoly,
}

impl Debug for DCRTPolyMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix")
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("params", &self.params)
            .field("entries", &self.entries.len())
            .finish()
    }
}

impl PartialEq for DCRTPolyMatrix {
    fn eq(&self, other: &Self) -> bool {
        if self.nrow != other.nrow || self.ncol != other.ncol {
            return false;
        }
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                let self_val = self.entry(i, j);
                let other_val = other.entry(i, j);
                if self_val != other_val {
                    return false;
                }
            }
        }
        true
    }
}

impl Eq for DCRTPolyMatrix {}

// Add getter methods for inner and params
impl DCRTPolyMatrix {
    pub fn entries(&self) -> &HashMap<(usize, usize), DCRTPoly> {
        &self.entries
    }

    pub fn params(&self) -> &DCRTPolyParams {
        &self.params
    }

    pub fn from_dense(
        params: &DCRTPolyParams,
        inner: Vec<Vec<DCRTPoly>>,
        nrow: usize,
        ncol: usize,
    ) -> Self {
        let zero_elem = DCRTPoly::const_zero(params);
        let mut entries = HashMap::new();

        for i in 0..nrow {
            for j in 0..ncol {
                if inner[i][j] != zero_elem {
                    entries.insert((i, j), inner[i][j].clone());
                }
            }
        }

        Self { entries, params: params.clone(), nrow, ncol, zero_elem }
    }

    pub fn to_dense(&self) -> Vec<Vec<DCRTPoly>> {
        let mut result = Vec::with_capacity(self.nrow);
        for i in 0..self.nrow {
            let mut row = Vec::with_capacity(self.ncol);
            for j in 0..self.ncol {
                if let Some(val) = self.entries.get(&(i, j)) {
                    row.push(val.clone());
                } else {
                    row.push(self.zero_elem.clone());
                }
            }
            result.push(row);
        }
        result
    }

    pub fn set_entry(&mut self, i: usize, j: usize, value: DCRTPoly) {
        if value == self.zero_elem {
            self.entries.remove(&(i, j));
        } else {
            self.entries.insert((i, j), value);
        }
    }
}

impl PolyMatrix for DCRTPolyMatrix {
    type P = DCRTPoly;

    fn from_poly_vec(params: &DCRTPolyParams, vec: Vec<Vec<DCRTPoly>>) -> Self {
        let nrow = vec.len();
        let ncol = vec[0].len();
        DCRTPolyMatrix::from_dense(params, vec, nrow, ncol)
    }

    fn entry(&self, i: usize, j: usize) -> &Self::P {
        self.entries.get(&(i, j)).unwrap_or(&self.zero_elem)
    }

    fn get_row(&self, i: usize) -> Vec<Self::P> {
        let mut row = Vec::with_capacity(self.ncol);
        for j in 0..self.ncol {
            row.push(self.entry(i, j).clone());
        }
        row
    }

    fn get_column(&self, j: usize) -> Vec<DCRTPoly> {
        let mut col = Vec::with_capacity(self.nrow);
        for i in 0..self.nrow {
            col.push(self.entry(i, j).clone());
        }
        col
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
        let mut entries = HashMap::new();
        for i in row_start..row_end {
            for j in column_start..column_end {
                if let Some(val) = self.entries.get(&(i, j)) {
                    entries.insert((i - row_start, j - column_start), val.clone());
                }
            }
        }
        Self { entries, params: self.params.clone(), nrow, ncol, zero_elem: self.zero_elem.clone() }
    }

    fn zero(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        Self {
            entries: HashMap::new(),
            params: params.clone(),
            nrow,
            ncol,
            zero_elem: DCRTPoly::const_zero(params),
        }
    }

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        let scalar = scalar.unwrap_or_else(|| DCRTPoly::const_one(params));
        let zero_elem = DCRTPoly::const_zero(params);
        let mut entries = HashMap::new();
        if scalar != zero_elem {
            for i in 0..size {
                entries.insert((i, i), scalar.clone());
            }
        }

        Self { entries, params: params.clone(), nrow: size, ncol: size, zero_elem }
    }

    fn transpose(&self) -> Self {
        let mut entries = HashMap::new();
        for ((i, j), val) in &self.entries {
            entries.insert((*j, *i), val.clone());
        }

        Self {
            entries,
            params: self.params.clone(),
            nrow: self.ncol,
            ncol: self.nrow,
            zero_elem: self.zero_elem.clone(),
        }
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

        let mut column_offsets = vec![0usize];
        let mut total_cols = self.ncol;

        for other in others {
            column_offsets.push(total_cols);
            total_cols += other.ncol;
        }

        let mut entries = HashMap::new();

        // Add entries from self
        for ((i, j), val) in &self.entries {
            entries.insert((*i, *j), val.clone());
        }

        // Add entries from others with column offset
        for (idx, other) in others.iter().enumerate() {
            let col_offset = column_offsets[idx + 1];
            for ((i, j), val) in &other.entries {
                entries.insert((*i, j + col_offset), val.clone());
            }
        }

        Self {
            entries,
            params: self.params.clone(),
            nrow: self.nrow,
            ncol: total_cols,
            zero_elem: self.zero_elem.clone(),
        }
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

        let mut row_offsets = vec![0usize];
        let mut total_rows = self.nrow;

        for other in others {
            row_offsets.push(total_rows);
            total_rows += other.nrow;
        }

        let mut entries = HashMap::new();

        // Add entries from self
        for ((i, j), val) in &self.entries {
            entries.insert((*i, *j), val.clone());
        }

        // Add entries from others with row offset
        for (idx, other) in others.iter().enumerate() {
            let row_offset = row_offsets[idx + 1];
            for ((i, j), val) in &other.entries {
                entries.insert((i + row_offset, *j), val.clone());
            }
        }

        Self {
            entries,
            params: self.params.clone(),
            nrow: total_rows,
            ncol: self.ncol,
            zero_elem: self.zero_elem.clone(),
        }
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[&Self]) -> Self {
        let mut row_offsets = vec![0usize];
        let mut col_offsets = vec![0usize];
        let mut total_rows = self.nrow;
        let mut total_cols = self.ncol;

        for other in others {
            row_offsets.push(total_rows);
            col_offsets.push(total_cols);
            total_rows += other.nrow;
            total_cols += other.ncol;
        }

        let mut entries = HashMap::new();

        // Add entries from self
        for ((i, j), val) in &self.entries {
            entries.insert((*i, *j), val.clone());
        }

        // Add entries from others with offsets
        for (idx, other) in others.iter().enumerate() {
            let row_offset = row_offsets[idx + 1];
            let col_offset = col_offsets[idx + 1];
            for ((i, j), val) in &other.entries {
                entries.insert((i + row_offset, j + col_offset), val.clone());
            }
        }

        Self {
            entries,
            params: self.params.clone(),
            nrow: total_rows,
            ncol: total_cols,
            zero_elem: self.zero_elem.clone(),
        }
    }

    fn tensor(&self, other: &Self) -> Self {
        let nrow = self.nrow * other.nrow;
        let ncol = self.ncol * other.ncol;
        let mut entries = HashMap::new();

        // For each entry in self
        for ((i1, j1), val1) in &self.entries {
            // For each entry in other
            for ((i2, j2), val2) in &other.entries {
                let new_i = i1 * other.nrow + i2;
                let new_j = j1 * other.ncol + j2;
                entries.insert((new_i, new_j), val1.clone() * val2.clone());
            }
        }

        Self { entries, params: self.params.clone(), nrow, ncol, zero_elem: self.zero_elem.clone() }
    }

    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self {
        let bit_length = params.modulus_bits();

        // Create gadget vector [1, 2, 2^2, ..., 2^(bit_length-1)]
        let mut poly_vec = Vec::with_capacity(bit_length);
        for i in 0u32..(bit_length as u32) {
            let value = BigInt::from(2).pow(i);
            poly_vec.push(DCRTPoly::from_const(params, &FinRingElem::new(value, params.modulus())));
        }
        let gadget_vector = Self::from_poly_vec(params, vec![poly_vec]);
        let identity = Self::identity(params, size, None);

        // Tensor product of identity and gadget vector
        identity.tensor(&gadget_vector)
    }

    fn decompose(&self) -> Self {
        let bit_length = self.params.modulus_bits();

        let mut decomposed_entries = HashMap::new();
        let decomposed_nrow = self.nrow * bit_length;

        info!("decompose decomposed_nrow{}", decomposed_nrow);

        for i in 0..self.nrow {
            for j in 0..self.ncol {
                let entry = self.entries.get(&(i, j)).unwrap_or(&self.zero_elem);
                let decomposition = entry.decompose(&self.params);

                for bit in 0..bit_length {
                    decomposed_entries
                        .insert((i * bit_length + bit, j), decomposition[bit].clone());
                }
            }
        }

        Self {
            entries: decomposed_entries,
            params: self.params.clone(),
            nrow: decomposed_nrow,
            ncol: self.ncol,
            zero_elem: self.zero_elem.clone(),
        }
    }

    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self {
        let mut new_entries = HashMap::new();
        for (&(i, j), entry) in &self.entries {
            new_entries.insert((i, j), entry.modulus_switch(&self.params, new_modulus.clone()));
        }

        Self {
            entries: new_entries,
            params: self.params.clone(),
            nrow: self.nrow,
            ncol: self.ncol,
            zero_elem: self.zero_elem.clone(),
        }
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

    // todo: if sparse this no need
    fn from_compact_bytes(
        params: &<Self::P as Poly>::Params,
        bytes: Vec<Bytes>,
        offset: usize,
    ) -> Self {
        let metadata = &bytes[0];

        let ring_dimension =
            u32::from_le_bytes([metadata[0], metadata[1], metadata[2], metadata[3]]);
        let nrow =
            u32::from_le_bytes([metadata[4], metadata[5], metadata[6], metadata[7]]) as usize;
        let ncol =
            u32::from_le_bytes([metadata[8], metadata[9], metadata[10], metadata[11]]) as usize;
        let byte_size =
            u32::from_le_bytes([metadata[12], metadata[13], metadata[14], metadata[15]]) as usize;

        assert_eq!(
            ring_dimension,
            params.ring_dimension(),
            "Ring dimension mismatch: {} != {}",
            ring_dimension,
            params.ring_dimension()
        );

        let modulus_bytes = params.modulus().to_bytes_le();
        assert!(
            byte_size <= modulus_bytes.len(),
            "byte_size must not be greater than modulus_bytes: {} > {}",
            byte_size,
            modulus_bytes.len()
        );

        let mut result = Self::zero(params, nrow, ncol);

        let mut idx = 1; // Start from 1 because 0 is metadata

        // for i in 0..nrow {
        //     for j in 0..ncol {
        //         let poly_bytes = &bytes[idx];
        //         idx += 1;
        //         let poly = DCRTPoly::from_compact_bytes(params, poly_bytes, offset);
        //         result.inner[i][j] = poly;
        //     }
        // }
        result
    }

    // todo: if sparse this no need
    fn to_compact_bytes(&self, byte_size: usize, offset: usize) -> Vec<Bytes> {
        let modulus_bytes = self.params.modulus().to_bytes_le();

        assert!(
            byte_size <= modulus_bytes.len(),
            "byte_size must not be greater than modulus_bytes: {} > {}",
            byte_size,
            modulus_bytes.len()
        );

        let mut result = Vec::new();

        let ring_dimension: u32 = self.params.ring_dimension();
        let nrow = self.nrow;
        let ncol = self.ncol;

        let mut metadata = Vec::new();
        metadata.extend_from_slice(&(ring_dimension).to_le_bytes());
        metadata.extend_from_slice(&(nrow as u32).to_le_bytes());
        metadata.extend_from_slice(&(ncol as u32).to_le_bytes());
        metadata.extend_from_slice(&(byte_size as u32).to_le_bytes());
        metadata.extend_from_slice(&(offset as u32).to_le_bytes());
        result.push(Bytes::from(metadata));

        // for i in 0..self.nrow {
        //     for j in 0..self.ncol {
        //         let poly = &self.inner[i][j];
        //         result.push(poly.to_compact_bytes(byte_size, offset));
        //     }
        // }
        result
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

    fn add(mut self, rhs: &DCRTPolyMatrix) -> Self::Output {
        #[cfg(debug_assertions)]
        if self.nrow != rhs.nrow || self.ncol != rhs.ncol {
            panic!(
                "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
                self.nrow, self.ncol, rhs.nrow, rhs.ncol
            );
        }

        for (&(i, j), rhs_value) in &rhs.entries {
            if let Some(self_value) = self.entries.get_mut(&(i, j)) {
                *self_value += rhs_value.clone();
            } else {
                self.entries.insert((i, j), rhs_value.clone());
            }
        }

        self
    }
}

impl Neg for DCRTPolyMatrix {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut negated_entries = HashMap::new();
        for ((i, j), value) in self.entries {
            negated_entries.insert((i, j), -value);
        }

        DCRTPolyMatrix {
            entries: negated_entries,
            params: self.params,
            nrow: self.nrow,
            ncol: self.ncol,
            zero_elem: self.zero_elem,
        }
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

        let mut result_entries = HashMap::new();
        let zero_elem = self.zero_elem.clone();

        for i in 0..nrow {
            for j in 0..ncol {
                let mut sum = zero_elem.clone();
                let mut has_contribution = false;
                for k in 0..self.ncol {
                    let self_value = self.entries.get(&(i, k)).unwrap_or(&self.zero_elem);
                    let rhs_value = rhs.entries.get(&(k, j)).unwrap_or(&rhs.zero_elem);
                    if !(self_value == &zero_elem) && !(rhs_value == &zero_elem) {
                        sum += self_value * rhs_value;
                        has_contribution = true;
                    }
                }
                if has_contribution && !(sum == zero_elem) {
                    result_entries.insert((i, j), sum);
                }
            }
        }

        DCRTPolyMatrix {
            entries: result_entries,
            params: self.params.clone(),
            nrow,
            ncol,
            zero_elem: self.zero_elem.clone(),
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
        let mut result_entries = HashMap::new();
        for (&(i, j), value) in &self.entries {
            let product = value * rhs;
            if !(product == self.zero_elem) {
                result_entries.insert((i, j), product);
            }
        }

        DCRTPolyMatrix {
            entries: result_entries,
            params: self.params.clone(),
            nrow: self.nrow,
            ncol: self.ncol,
            zero_elem: self.zero_elem.clone(),
        }
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

    fn sub(mut self, rhs: &DCRTPolyMatrix) -> Self::Output {
        #[cfg(debug_assertions)]
        if self.nrow != rhs.nrow || self.ncol != rhs.ncol {
            panic!(
                "Subtraction requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
                self.nrow, self.ncol, rhs.nrow, rhs.ncol
            );
        }
        for (&(i, j), rhs_value) in &rhs.entries {
            if let Some(self_value) = self.entries.get_mut(&(i, j)) {
                *self_value -= rhs_value.clone();
                if self_value == &self.zero_elem {
                    self.entries.remove(&(i, j));
                }
            } else {
                self.entries.insert((i, j), -rhs_value.clone());
            }
        }

        self
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use num_bigint::BigUint;

    use super::*;
    use crate::poly::{
        dcrt::{DCRTPolyParams, DCRTPolyUniformSampler},
        sampler::PolyUniformSampler,
    };

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
        matrix.set_entry(0, 0, DCRTPoly::from_const(&params, &value));
        matrix.set_entry(1, 1, DCRTPoly::from_const(&params, &value));
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
        matrix1.set_entry(0, 0, DCRTPoly::from_const(&params, &value));
        matrix1.set_entry(1, 1, DCRTPoly::from_const(&params, &value));

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
        matrix1.set_entry(0, 0, DCRTPoly::from_const(&params, &value));

        let mut matrix2 = DCRTPolyMatrix::zero(&params, 2, 2);
        matrix2.set_entry(1, 1, DCRTPoly::from_const(&params, &value));

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
        matrix1.set_entry(0, 0, DCRTPoly::from_const(&params, &value));

        let mut matrix2 = DCRTPolyMatrix::zero(&params, 2, 2);
        matrix2.set_entry(0, 0, DCRTPoly::from_const(&params, &value));

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

    #[test]
    fn test_mul_tensor_identity() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix S (2x12)
        let s = sampler.sample_uniform(&params, 2, 12, crate::poly::sampler::DistType::FinRingDist);

        // Create 'other' matrix (3x3)
        let other =
            sampler.sample_uniform(&params, 3, 3, crate::poly::sampler::DistType::FinRingDist);

        // Perform S * (I_4 ⊗ other)
        let result = s.mul_tensor_identity(&other, 4);

        // Check dimensions
        assert_eq!(result.row_size(), 2);
        assert_eq!(result.col_size(), 12);

        let identity = DCRTPolyMatrix::identity(&params, 4, None);

        // Check result
        let expected_result = s * (identity.tensor(&other));

        assert_eq!(expected_result.row_size(), 2);
        assert_eq!(expected_result.col_size(), 12);
        assert_eq!(result, expected_result)
    }

    #[test]
    #[should_panic(expected = "value_bytes exceeds the specified byte_size:")]
    fn test_to_compact_bytes_failure_1() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix (2x12)
        let mat =
            sampler.sample_uniform(&params, 2, 12, crate::poly::sampler::DistType::FinRingDist);

        // Choose a byte size that is less than the bytes necessary to represent the coefficients
        let byte_size = 1;
        mat.to_compact_bytes(byte_size, 0);
    }

    #[test]
    #[should_panic(expected = "byte_size must not be greater than modulus_bytes")]
    fn test_to_compact_bytes_failure_2() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix (2x12)
        let mat =
            sampler.sample_uniform(&params, 2, 12, crate::poly::sampler::DistType::FinRingDist);

        // Choose a bit size that is greater than the modulus bit size
        let byte_size = params.modulus().to_bytes_le().len() + 1;
        mat.to_compact_bytes(byte_size, 0);
    }

    #[test]
    fn test_to_compact_bytes() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        let nrow = 2;
        let ncol = 12;
        let ring_dimension = params.ring_dimension();

        // Create matrix (2x12)
        let mat =
            sampler.sample_uniform(&params, nrow, ncol, crate::poly::sampler::DistType::BitDist);

        let byte_size = 1;
        let bytes = mat.to_compact_bytes(byte_size, 0);

        // the vector should contain 1 (metadata) + nrow * ncol elements
        assert_eq!(bytes.len(), 1 + (nrow * ncol));

        // the total byte size should be 4 * 5 (metadata) + nrow * ncol * ring_dimension * byte_size
        // (coefficients)
        let total_byte_size: usize = bytes.iter().map(|b| b.len()).sum();
        assert_eq!(total_byte_size, (4 * 5) + (nrow * ncol * ring_dimension as usize * byte_size));
    }

    #[test]
    fn test_from_compact_bytes() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        let nrow = 2;
        let ncol = 12;

        let mat = sampler.sample_uniform(
            &params,
            nrow,
            ncol,
            crate::poly::sampler::DistType::FinRingDist,
        );

        let byte_size = params.modulus().to_bytes_le().len();

        let bytes = mat.to_compact_bytes(byte_size, 0);
        let new_mat = DCRTPolyMatrix::from_compact_bytes(&params, bytes, 0);

        assert_eq!(mat, new_mat);

        let bin_mat =
            sampler.sample_uniform(&params, nrow, ncol, crate::poly::sampler::DistType::BitDist);

        let byte_size = 1;

        let bytes = bin_mat.to_compact_bytes(byte_size, 0);
        let new_bin_mat = DCRTPolyMatrix::from_compact_bytes(&params, bytes, 0);

        assert_eq!(bin_mat, new_bin_mat);

        let sigma = 4.57825;

        let gauss_mat = sampler.sample_uniform(
            &params,
            nrow,
            ncol,
            crate::poly::sampler::DistType::GaussDist { sigma },
        );

        let bound = 3.0 * sigma;

        let bound_ceil = bound.ceil() as u32;

        // byte_size is bytes necessary to represent bound*2 which is the maximum possible sampled
        // coefficient value
        let byte_size = (2 * bound_ceil).to_le_bytes().len();
        let bytes = gauss_mat.to_compact_bytes(byte_size, bound_ceil as usize);
        let new_gauss_mat = DCRTPolyMatrix::from_compact_bytes(&params, bytes, bound_ceil as usize);
        assert_eq!(gauss_mat, new_gauss_mat);
    }
}
