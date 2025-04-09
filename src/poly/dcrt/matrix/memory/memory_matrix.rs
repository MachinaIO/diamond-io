use crate::{parallel_iter, poly::MatrixElem};
use rayon::prelude::*;
use std::{
    env,
    fmt::Debug,
    ops::{Add, Mul, Neg, Range, Sub},
};

#[derive(Clone)]
pub struct MemoryMatrix<T: MatrixElem> {
    pub inner: Vec<Vec<T>>,
    pub params: T::Params,
    pub nrow: usize,
    pub ncol: usize,
}

impl<T: MatrixElem> Debug for MemoryMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix")
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("params", &self.params)
            .field("inner", &self.inner)
            .finish()
    }
}

impl<T: MatrixElem> PartialEq for MemoryMatrix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner && self.nrow == other.nrow && self.ncol == other.ncol
    }
}

impl<T: MatrixElem> Eq for MemoryMatrix<T> {}

// Add getter methods for inner and params
impl<T: MatrixElem> MemoryMatrix<T> {
    pub fn new_empty(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        let inner = vec![vec![T::new_empty(); ncol]; nrow];
        Self { inner, params: params.clone(), nrow, ncol }
    }

    pub fn inner(&self) -> &Vec<Vec<T>> {
        &self.inner
    }

    pub fn replace_entries<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
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

    pub fn entry(&self, i: usize, j: usize) -> T {
        self.inner[i][j].clone()
    }

    pub fn get_row(&self, i: usize) -> Vec<T> {
        self.inner[i].clone()
    }

    pub fn get_column(&self, j: usize) -> Vec<T> {
        self.inner.iter().map(|row| row[j].clone()).collect()
    }

    pub fn size(&self) -> (usize, usize) {
        (self.nrow, self.ncol)
    }

    pub fn row_size(&self) -> usize {
        self.nrow
    }

    pub fn col_size(&self) -> usize {
        self.ncol
    }

    pub fn slice(
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

        Self { inner: c, params: self.params.clone(), nrow, ncol }
    }

    pub fn block_entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        rows.into_par_iter()
            .map(|i| {
                // Use slice indexing on the row. The slice is converted into a Vec by `.to_vec()`.
                self.inner[i][cols.clone()].to_vec()
            })
            .collect()
    }

    pub fn zero(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        let mut c = Vec::with_capacity(nrow);
        let zero_elem = T::zero(params);
        for _ in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for _ in 0..ncol {
                row.push(zero_elem.clone());
            }
            c.push(row);
        }
        Self { inner: c, params: params.clone(), nrow, ncol }
    }

    pub fn identity(params: &T::Params, size: usize, scalar: Option<T>) -> Self {
        let nrow = size;
        let ncol = size;
        let scalar = scalar.unwrap_or_else(|| T::one(params));
        let zero_elem = T::zero(params);
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

        Self { inner: result, params: params.clone(), nrow, ncol }
    }

    pub fn transpose(&self) -> Self {
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

        Self { inner: result, params: self.params.clone(), nrow, ncol }
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    pub fn concat_columns(&self, others: &[&Self]) -> Self {
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
        let result: Vec<Vec<T>> = parallel_iter!(0..self.nrow)
            .map(|i| {
                let mut row: Vec<&T> = Vec::with_capacity(ncol);
                row.extend(self.inner[i].iter());
                for other in others {
                    row.extend(other.inner[i].iter());
                }
                row.into_iter().cloned().collect()
            })
            .collect();

        Self { inner: result, params: self.params.clone(), nrow: self.nrow, ncol }
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    pub fn concat_rows(&self, others: &[&Self]) -> Self {
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

        Self { inner: result, params: self.params.clone(), nrow, ncol: self.ncol }
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    pub fn concat_diag(&self, others: &[&Self]) -> Self {
        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();

        let zero_elem = T::zero(&self.params);

        let mut result: Vec<Vec<T>> = Vec::with_capacity(nrow);

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
                    .collect::<Vec<Vec<T>>>()
                    .into_iter(),
            );
            col_offset += other.ncol;
        }

        Self { inner: result, params: self.params.clone(), nrow, ncol }
    }

    pub fn tensor(&self, other: &Self) -> Self {
        let nrow = self.nrow * other.nrow;
        let ncol = self.ncol * other.ncol;
        let mut result: Vec<Vec<T>> = Vec::with_capacity(nrow);
        for i1 in 0..self.nrow {
            for i2 in 0..other.nrow {
                let mut row = Vec::with_capacity(ncol);
                for j1 in 0..self.ncol {
                    for j2 in 0..other.ncol {
                        row.push(self.inner[i1][j1].clone() * &other.inner[i2][j2]);
                    }
                }
                result.push(row);
            }
        }

        Self { inner: result, params: self.params.clone(), nrow, ncol }
    }
}

// ====== Arithmetic ======

impl<T: MatrixElem> Add for MemoryMatrix<T> {
    type Output = MemoryMatrix<T>;
    fn add(self, rhs: Self) -> Self {
        self + &rhs
    }
}

// Implement addition of a matrix by a matrix reference
impl<T: MatrixElem> Add<&MemoryMatrix<T>> for MemoryMatrix<T> {
    type Output = Self;

    fn add(self, rhs: &MemoryMatrix<T>) -> Self::Output {
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

impl<T: MatrixElem> Neg for MemoryMatrix<T> {
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

        Self { inner: c, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}

impl<T: MatrixElem> Mul for MemoryMatrix<T> {
    type Output = MemoryMatrix<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&MemoryMatrix<T>> for MemoryMatrix<T> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<T: MatrixElem> Mul<&MemoryMatrix<T>> for &MemoryMatrix<T> {
    type Output = MemoryMatrix<T>;

    fn mul(self, rhs: &MemoryMatrix<T>) -> Self::Output {
        let nrow = self.nrow;
        let ncol = rhs.ncol;
        #[cfg(debug_assertions)]
        if rhs.nrow != self.ncol {
            panic!(
                "Multiplication condition failed: rhs.nrow ({}) must equal self.ncol ({})",
                rhs.nrow, self.ncol
            );
        }

        MemoryMatrix {
            inner: parallel_iter!(0..nrow)
                .map(|i| {
                    parallel_iter!(0..ncol)
                        .map(|j| {
                            let mut sum = self.inner[i][0].clone() * &rhs.inner[0][j];
                            for k in 1..self.ncol {
                                sum += self.inner[i][k].clone() * &rhs.inner[k][j];
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
impl<T: MatrixElem> Mul<T> for MemoryMatrix<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&T> for MemoryMatrix<T> {
    type Output = Self;

    fn mul(self, rhs: &T) -> Self::Output {
        &self * rhs
    }
}

impl<T: MatrixElem> Mul<&T> for &MemoryMatrix<T> {
    type Output = MemoryMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;

        let result: Vec<Vec<T>> = parallel_iter!(0..nrow)
            .map(|i| parallel_iter!(0..ncol).map(|j| self.inner[i][j].clone() * rhs).collect())
            .collect();

        MemoryMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }
}

// Implement subtraction for matrices
impl<T: MatrixElem> Sub for MemoryMatrix<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<T: MatrixElem> Sub<&MemoryMatrix<T>> for MemoryMatrix<T> {
    type Output = Self;

    fn sub(self, rhs: &MemoryMatrix<T>) -> Self::Output {
        &self - rhs
    }
}

// Implement subtraction of a matrix by a matrix reference
impl<T: MatrixElem> Sub<&MemoryMatrix<T>> for &MemoryMatrix<T> {
    type Output = MemoryMatrix<T>;

    fn sub(self, rhs: &MemoryMatrix<T>) -> Self::Output {
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

        MemoryMatrix { inner: result, params: self.params.clone(), ncol, nrow }
    }
}
