use crate::{
    parallel_iter,
    poly::{MatrixElem, MatrixParams},
    utils::block_size,
};
use itertools::Itertools;
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Range, Sub},
};

#[derive(Clone)]
pub struct BaseMatrix<T: MatrixElem> {
    pub params: T::Params,
    pub inner: Vec<Vec<T>>,
    pub nrow: usize,
    pub ncol: usize,
}

impl<T: MatrixElem> BaseMatrix<T> {
    pub fn new_empty(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        let inner = vec![vec![T::zero(params); ncol]; nrow];
        Self { params: params.clone(), inner, nrow, ncol }
    }

    pub fn entry_size(&self) -> usize {
        self.params.entry_size()
    }

    pub fn block_entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        parallel_iter!(rows)
            .map(|i| self.inner[i][cols.start..cols.end].to_vec())
            .collect::<Vec<Vec<_>>>()
    }

    pub fn replace_entries<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        if self.nrow == 0 || self.ncol == 0 {
            return;
        }
        let polys = f(rows.clone(), cols.clone());
        debug_assert_eq!(polys.len(), rows.len());
        debug_assert_eq!(polys[0].len(), cols.len());
        self.inner[rows.start..rows.end].par_iter_mut().enumerate().for_each(|(i, row_data)| {
            row_data[cols.start..cols.end].clone_from_slice(&polys[i]);
        });
    }

    pub fn replace_entries_diag<F>(&mut self, diags: Range<usize>, f: F)
    where
        F: Fn(Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let polys = f(diags.clone());
        debug_assert_eq!(polys.len(), diags.len());
        debug_assert_eq!(polys[0].len(), diags.len());
        self.inner[diags.start..diags.end].par_iter_mut().enumerate().for_each(|(i, row_data)| {
            row_data[diags.start..diags.end].clone_from_slice(&polys[i]);
        });
    }

    pub fn replace_entries_with_expand<F>(
        &mut self,
        rows: Range<usize>,
        cols: Range<usize>,
        row_scale: usize,
        col_scale: usize,
        f: F,
    ) where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let polys = f(rows.clone(), cols.clone());
        debug_assert_eq!(polys.len(), rows.len() * row_scale);
        debug_assert_eq!(polys[0].len(), cols.len() * col_scale);
        self.inner[rows.start * row_scale..rows.end * row_scale]
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, row_data)| {
                row_data[cols.start * col_scale..cols.end * col_scale].clone_from_slice(&polys[i]);
            });
    }

    pub fn entry(&self, i: usize, j: usize) -> T {
        self.block_entries(i..i + 1, j..j + 1)[0][0].clone()
    }

    pub fn get_row(&self, i: usize) -> Vec<T> {
        self.block_entries(i..i + 1, 0..self.ncol)[0].clone()
    }

    pub fn get_column(&self, j: usize) -> Vec<T> {
        self.block_entries(0..self.nrow, j..j + 1).iter().map(|row| row[0].clone()).collect()
    }

    pub fn size(&self) -> (usize, usize) {
        (self.nrow, self.ncol)
    }

    pub fn slice(
        &self,
        row_start: usize,
        row_end: usize,
        col_start: usize,
        col_end: usize,
    ) -> Self {
        let nrow = row_end - row_start;
        let ncol = col_end - col_start;

        let inner: Vec<_> = (row_start..row_end)
            .into_par_iter()
            .map(|i| {
                (col_start..col_end)
                    .into_par_iter()
                    .map(|j| self.inner[i][j].clone())
                    .collect::<Vec<_>>()
            })
            .collect();

        Self { inner, params: self.params.clone(), nrow, ncol }
    }

    pub fn zero(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        Self::new_empty(params, nrow, ncol)
    }

    pub fn identity(params: &T::Params, size: usize, scalar: Option<T>) -> Self {
        let nrow = size;
        let ncol = size;
        let scalar = scalar.unwrap_or_else(|| T::one(params));
        let zero_elem = T::zero(params);
        let inner: Vec<Vec<T>> = (0..size)
            .into_par_iter() //
            .map(|i| {
                (0..size)
                    .into_par_iter()
                    .map(|j| if i == j { scalar.clone() } else { zero_elem.clone() })
                    .collect()
            })
            .collect();

        Self { inner, params: params.clone(), nrow, ncol }
    }

    pub fn transpose(&self) -> Self {
        let nrow = self.ncol;
        let ncol = self.nrow;
        let inner: Vec<Vec<T>> = (0..self.ncol)
            .into_par_iter()
            .map(|i| {
                (0..self.nrow).into_par_iter().map(|j| self.inner[j][i].clone()).collect::<Vec<T>>()
            })
            .collect();

        Self { inner, params: self.params.clone(), nrow, ncol }
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
                    "Concat error: while the shape of the first matrix is ({}, {}), that of the {}-th matrix is ({},{})",
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
        let nrow_total = self.nrow * other.nrow;
        let ncol_total = self.ncol * other.ncol;

        let index_pairs: Vec<(usize, usize)> =
            (0..self.nrow).flat_map(|i1| (0..other.nrow).map(move |i2| (i1, i2))).collect();

        let result: Vec<Vec<T>> = index_pairs
            .into_par_iter()
            .map(|(i1, i2)| {
                let mut row = Vec::with_capacity(ncol_total);
                for j1 in 0..self.ncol {
                    for j2 in 0..other.ncol {
                        row.push(self.inner[i1][j1].clone() * &other.inner[i2][j2]);
                    }
                }
                row
            })
            .collect();

        Self { params: self.params.clone(), inner: result, nrow: nrow_total, ncol: ncol_total }
    }

    #[inline]
    pub fn set_entry(&mut self, i: usize, j: usize, elem: T) {
        debug_assert!(i < self.nrow, "row index {i} ≥ nrow = {}", self.nrow);
        debug_assert!(j < self.ncol, "col index {j} ≥ ncol = {}", self.ncol);

        self.inner[i][j] = elem;
    }
}

impl<T: MatrixElem> Debug for BaseMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let fmt = f
            .debug_struct("BaseMatrix")
            .field("params", &self.params)
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("inner", &self.inner)
            .finish();
        fmt
    }
}

impl<T: MatrixElem> PartialEq for BaseMatrix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner &&
            self.params == other.params &&
            self.nrow == other.nrow &&
            self.ncol == other.ncol
    }
}

impl<T: MatrixElem> Eq for BaseMatrix<T> {}

impl<T: MatrixElem> Add for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<T: MatrixElem> Add<&BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: &Self) -> Self::Output {
        &self + rhs
    }
}

impl<T: MatrixElem> Add<&BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: &BaseMatrix<T>) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );

        let mut new_matrix = BaseMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            add_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MatrixElem> Sub<BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<T: MatrixElem> Sub<&BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: &Self) -> Self::Output {
        &self - rhs
    }
}

impl<T: MatrixElem> Sub<&BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: &BaseMatrix<T>) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );

        let mut new_matrix = BaseMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            sub_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MatrixElem> Mul<BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<T: MatrixElem> Mul<BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: BaseMatrix<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &BaseMatrix<T>) -> Self::Output {
        debug_assert!(
            self.ncol == rhs.nrow,
            "Multiplication condition failed: self.ncol ({}) must equal rhs.nrow ({})",
            self.ncol,
            rhs.nrow
        );

        let mut new_matrix = BaseMatrix::new_empty(&self.params, self.nrow, rhs.ncol);
        let (_, ip_offsets) = block_offsets(0..0, 0..self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            ip_offsets
                .iter()
                .tuple_windows()
                .map(|(cur_block_ip_idx, next_block_ip_idx)| {
                    let self_block_polys = self
                        .block_entries(row_offsets.clone(), *cur_block_ip_idx..*next_block_ip_idx);
                    let other_block_polys = rhs
                        .block_entries(*cur_block_ip_idx..*next_block_ip_idx, col_offsets.clone());
                    mul_block_matrices(self_block_polys, other_block_polys)
                })
                .reduce(|acc, muled| add_block_matrices(muled, &acc))
                .unwrap()
        };
        new_matrix.replace_entries(0..self.nrow, 0..rhs.ncol, f);
        new_matrix
    }
}

impl<T: MatrixElem> Mul<T> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&T> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        &self * rhs
    }
}

impl<T: MatrixElem> Mul<T> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&T> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        let mut new_matrix = BaseMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
                .into_iter()
                .map(|row| row.into_iter().map(|elem| elem * rhs).collect::<Vec<T>>())
                .collect::<Vec<Vec<T>>>()
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MatrixElem> Neg for BaseMatrix<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut new_matrix = BaseMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
                .into_iter()
                .map(|row| row.into_iter().map(|elem| -elem.clone()).collect())
                .collect()
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

pub fn block_offsets(rows: Range<usize>, cols: Range<usize>) -> (Vec<usize>, Vec<usize>) {
    let block_size = block_size();
    block_offsets_distinct_block_sizes(rows, cols, block_size, block_size)
}

pub fn block_offsets_distinct_block_sizes(
    rows: Range<usize>,
    cols: Range<usize>,
    row_block_size: usize,
    col_block_size: usize,
) -> (Vec<usize>, Vec<usize>) {
    let nrow = rows.end - rows.start;
    let ncol = cols.end - cols.start;
    let num_blocks_row = nrow.div_ceil(row_block_size);
    let num_blocks_col = ncol.div_ceil(col_block_size);
    let mut row_offsets = vec![rows.start];
    for _ in 0..num_blocks_row {
        let last_row_offset = row_offsets.last().unwrap();
        let sub = rows.end - last_row_offset;
        let len = if sub < row_block_size { sub } else { row_block_size };
        row_offsets.push(last_row_offset + len);
    }
    let mut col_offsets = vec![cols.start];
    for _ in 0..num_blocks_col {
        let last_col_offset = col_offsets.last().unwrap();
        let sub = cols.end - last_col_offset;
        let len = if sub < col_block_size { sub } else { col_block_size };
        col_offsets.push(last_col_offset + len);
    }
    (row_offsets, col_offsets)
}

fn add_block_matrices<T: MatrixElem>(lhs: Vec<Vec<T>>, rhs: &[Vec<T>]) -> Vec<Vec<T>> {
    let nrow = lhs.len();
    let ncol = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol).map(|j| lhs[i][j].clone() + &rhs[i][j]).collect::<Vec<T>>()
        })
        .collect::<Vec<Vec<T>>>()
}

fn sub_block_matrices<T: MatrixElem>(lhs: Vec<Vec<T>>, rhs: &[Vec<T>]) -> Vec<Vec<T>> {
    let nrow = lhs.len();
    let ncol = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol).map(|j| lhs[i][j].clone() - &rhs[i][j]).collect::<Vec<T>>()
        })
        .collect::<Vec<Vec<T>>>()
}

fn mul_block_matrices<T: MatrixElem>(lhs: Vec<Vec<T>>, rhs: Vec<Vec<T>>) -> Vec<Vec<T>> {
    let nrow = lhs.len();
    let ncol = rhs[0].len();
    let n_inner = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol)
                .map(|j: usize| {
                    (0..n_inner)
                        .map(|k| lhs[i][k].clone() * &rhs[k][j])
                        .reduce(|acc, prod| acc + prod)
                        .unwrap()
                })
                .collect::<Vec<T>>()
        })
        .collect::<Vec<Vec<T>>>()
}
