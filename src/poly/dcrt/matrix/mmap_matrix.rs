use crate::parallel_iter;
use itertools::Itertools;
use memmap2::{Mmap, MmapMut, MmapOptions};
// use once_cell::sync::OnceCell;
use rayon::prelude::*;
use std::{
    env,
    fmt::Debug,
    fs::File,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Range, Sub, SubAssign},
};
// use sysinfo::System;
use tempfile::tempfile;

// static BLOCK_SIZE: OnceCell<usize> = OnceCell::new();

pub trait MmapMatrixParams: Debug + Clone + PartialEq + Eq + Send + Sync {
    fn entry_size(&self) -> usize;
}

pub trait MmapMatrixElem:
    Sized
    + Clone
    + Debug
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + Send
    + Sync
{
    type Params: MmapMatrixParams;
    fn zero(params: &Self::Params) -> Self;
    fn one(params: &Self::Params) -> Self;
    fn from_bytes_to_elem(params: &Self::Params, bytes: &[u8]) -> Self;
    fn as_elem_to_bytes(&self) -> Vec<u8>;
}

pub struct MmapMatrix<T: MmapMatrixElem> {
    pub params: T::Params,
    pub file: File,
    pub nrow: usize,
    pub ncol: usize,
}

impl<T: MmapMatrixElem> MmapMatrix<T> {
    pub fn new_empty(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        let entry_size = params.entry_size();
        let len = entry_size * nrow * ncol;
        let file = tempfile().expect("failed to open file");
        file.set_len(len as u64).expect("failed to set file length");
        Self { params: params.clone(), file, nrow, ncol }
    }

    pub fn entry_size(&self) -> usize {
        self.params.entry_size()
    }

    pub fn block_entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        let entry_size = self.entry_size();
        parallel_iter!(rows)
            .map(|i| {
                let offset = entry_size * (i * self.ncol + cols.start);
                let mmap = map_file(&self.file, offset, entry_size * cols.len());
                let row_vec = mmap.to_vec();
                let row_col_vec = row_vec
                    .chunks(entry_size)
                    .map(|entry| T::from_bytes_to_elem(&self.params, entry))
                    .collect_vec();
                drop(mmap);
                row_col_vec
            })
            .collect::<Vec<Vec<_>>>()
    }

    pub fn replace_entries<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let (row_offsets, col_offsets) = block_offsets(rows.clone(), cols.clone());
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let new_entries = f(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            self.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_entries,
                            );
                        }
                    },
                );
            },
        );
    }

    pub fn replace_entries_diag<F>(&mut self, diags: Range<usize>, f: F)
    where
        F: Fn(Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let (offsets, _) = block_offsets(diags.clone(), 0..0);
        parallel_iter!(offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_diag_idx, next_block_diag_idx)| {
                let new_entries = f(*cur_block_diag_idx..*next_block_diag_idx);
                // This is secure because the modified entries are not overlapped among
                // threads
                unsafe {
                    self.replace_block_entries(
                        *cur_block_diag_idx..*next_block_diag_idx,
                        *cur_block_diag_idx..*next_block_diag_idx,
                        new_entries,
                    );
                }
            },
        );
    }

    pub(crate) unsafe fn replace_block_entries(
        &self,
        rows: Range<usize>,
        cols: Range<usize>,
        new_entries: Vec<Vec<T>>,
    ) {
        debug_assert_eq!(new_entries.len(), rows.end - rows.start);
        debug_assert_eq!(new_entries[0].len(), cols.end - cols.start);
        let entry_size = self.entry_size();
        let row_start = rows.start;
        parallel_iter!(rows).for_each(|i| {
            let offset = entry_size * (i * self.ncol + cols.start);
            let mut mmap = unsafe { map_file_mut(&self.file, offset, entry_size * cols.len()) };
            let bytes = new_entries[i - row_start]
                .iter()
                .flat_map(|poly| poly.as_elem_to_bytes())
                .collect::<Vec<_>>();
            mmap.copy_from_slice(&bytes);
            drop(mmap);
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
        let mut new_matrix = Self::new_empty(&self.params, nrow, ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let row_offsets = row_start + row_offsets.start..row_start + row_offsets.end;
            let col_offsets = col_start + col_offsets.start..col_start + col_offsets.end;

            self.block_entries(row_offsets, col_offsets)
        };
        new_matrix.replace_entries(0..nrow, 0..ncol, f);
        new_matrix
    }

    pub fn zero(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        Self::new_empty(params, nrow, ncol)
    }

    pub fn identity(params: &T::Params, size: usize, scalar: Option<T>) -> Self {
        let mut new_matrix = Self::new_empty(params, size, size);
        let scalar = scalar.unwrap_or_else(|| T::one(params));
        let f = |offsets: Range<usize>| -> Vec<Vec<T>> {
            let mut new_entries = vec![vec![T::zero(params); offsets.len()]; offsets.len()];
            for i in offsets {
                new_entries[i][i] = scalar.clone();
            }
            new_entries
        };
        new_matrix.replace_entries_diag(0..size, f);
        new_matrix
    }

    pub fn transpose(&self) -> Self {
        let mut new_matrix = Self::new_empty(&self.params, self.ncol, self.nrow);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let mut new_entries =
                vec![vec![T::zero(&self.params); col_offsets.len()]; row_offsets.len()];
            let cur_entries = self.block_entries(col_offsets.clone(), row_offsets.clone());
            for i in row_offsets {
                for j in col_offsets.clone() {
                    new_entries[i][j] = cur_entries[j][i].clone();
                }
            }
            new_entries
        };
        new_matrix.replace_entries(0..self.ncol, 0..self.nrow, f);
        new_matrix
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
        let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);
        let mut new_matrix = Self::new_empty(&self.params, self.nrow, updated_ncol);
        let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, self_f);
        let mut col_acc = self.ncol;
        for other in others.iter() {
            let other_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                let col_offsets = col_offsets.start - col_acc..col_offsets.end - col_acc;
                other.block_entries(row_offsets, col_offsets)
            };
            new_matrix.replace_entries(0..self.nrow, col_acc..col_acc + other.ncol, other_f);
            col_acc += other.ncol;
        }
        debug_assert_eq!(col_acc, updated_ncol);
        new_matrix
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
        let updated_nrow = others.iter().fold(self.nrow, |acc, other| acc + other.nrow);

        let mut new_matrix = Self::new_empty(&self.params, updated_nrow, self.ncol);
        let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, self_f);

        let mut row_acc = self.nrow;
        for other in others.iter() {
            let other_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                let row_offsets = row_offsets.start - row_acc..row_offsets.end - row_acc;
                other.block_entries(row_offsets, col_offsets)
            };
            new_matrix.replace_entries(row_acc..row_acc + other.nrow, 0..self.ncol, other_f);
            row_acc += other.nrow;
        }
        debug_assert_eq!(row_acc, updated_nrow);
        new_matrix
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    pub fn concat_diag(&self, others: &[&Self]) -> Self {
        let updated_nrow = others.iter().fold(self.nrow, |acc, other| acc + other.nrow);
        let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);

        let mut new_matrix = Self::new_empty(&self.params, updated_nrow, updated_ncol);
        let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, self_f);
        let mut row_acc = self.nrow;
        let mut col_acc = self.ncol;
        for other in others.iter() {
            let other_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                let row_offsets = row_offsets.start - row_acc..row_offsets.end - row_acc;
                let col_offsets = col_offsets.start - col_acc..col_offsets.end - col_acc;
                other.block_entries(row_offsets, col_offsets)
            };
            new_matrix.replace_entries(
                row_acc..row_acc + other.nrow,
                col_acc..col_acc + other.ncol,
                other_f,
            );
            row_acc += other.nrow;
            col_acc += other.ncol;
        }
        debug_assert_eq!(row_acc, updated_nrow);
        debug_assert_eq!(col_acc, updated_ncol);
        new_matrix
    }

    pub fn tensor(&self, other: &Self) -> Self {
        let new_matrix =
            Self::new_empty(&self.params, self.nrow * other.nrow, self.ncol * other.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..other.nrow, 0..other.ncol);
        parallel_iter!(0..self.nrow).for_each(|i| {
            parallel_iter!(0..self.ncol).for_each(|j| {
                let scalar = self.entry(i, j);
                let sub_matrix = other.clone() * scalar;

                parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_row_idx, next_block_row_idx)| {
                        parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                            |(cur_block_col_idx, next_block_col_idx)| {
                                let sub_block_polys = sub_matrix.block_entries(
                                    *cur_block_row_idx..*next_block_row_idx,
                                    *cur_block_col_idx..*next_block_col_idx,
                                );
                                // This is secure because the modified entries are not overlapped
                                // among threads
                                unsafe {
                                    new_matrix.replace_block_entries(
                                        i * sub_matrix.nrow + *cur_block_row_idx..
                                            i * sub_matrix.nrow + *next_block_row_idx,
                                        j * sub_matrix.ncol + *cur_block_col_idx..
                                            j * sub_matrix.ncol + *next_block_col_idx,
                                        sub_block_polys,
                                    );
                                }
                            },
                        );
                    },
                );
            });
        });
        new_matrix
    }
}

impl<T: MmapMatrixElem> Debug for MmapMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MmapMatrix")
            .field("params", &self.params)
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("file", &self.file)
            .finish()
    }
}

impl<T: MmapMatrixElem> Clone for MmapMatrix<T> {
    fn clone(&self) -> Self {
        let mut new_matrix = Self::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> PartialEq for MmapMatrix<T> {
    fn eq(&self, other: &Self) -> bool {
        if self.params != other.params || self.nrow != other.nrow || self.ncol != other.ncol {
            return false;
        }
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).all(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).all(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let other_block_polys = other.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        self_block_polys
                            .iter()
                            .zip(other_block_polys.iter())
                            .all(|(self_poly, other_poly)| self_poly == other_poly)
                    },
                )
            },
        )
    }
}

impl<T: MmapMatrixElem> Eq for MmapMatrix<T> {}

impl<T: MmapMatrixElem> Drop for MmapMatrix<T> {
    fn drop(&mut self) {
        self.file.set_len(0).expect("failed to truncate file");
    }
}

impl<T: MmapMatrixElem> Add for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<T: MmapMatrixElem> Add<&MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn add(self, rhs: &MmapMatrix<T>) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            add_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Sub for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<T: MmapMatrixElem> Sub<&MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn sub(self, rhs: &MmapMatrix<T>) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            sub_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Mul for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: MmapMatrixElem> Mul<&MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: &Self) -> Self::Output {
        debug_assert!(
            self.ncol == rhs.nrow,
            "Multiplication condition failed: self.ncol ({}) must equal rhs.nrow ({})",
            self.ncol,
            rhs.nrow
        );
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, rhs.ncol);
        let (_, ip_offsets) = block_offsets(0..0, 0..self.ncol);

        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let mut ip_sum =
                vec![vec![T::zero(&self.params); col_offsets.len()]; row_offsets.len()];
            for (cur_block_ip_idx, next_block_ip_idx) in ip_offsets.iter().tuple_windows() {
                let self_block_polys =
                    self.block_entries(row_offsets.clone(), *cur_block_ip_idx..*next_block_ip_idx);
                let other_block_polys =
                    rhs.block_entries(*cur_block_ip_idx..*next_block_ip_idx, col_offsets.clone());
                let muled = mul_block_matrices(
                    &self.params,
                    self_block_polys.clone(),
                    other_block_polys.clone(),
                );
                ip_sum = add_block_matrices(muled, &ip_sum);
            }
            ip_sum
        };
        new_matrix.replace_entries(0..self.nrow, 0..rhs.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Mul<T> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MmapMatrixElem> Mul<&T> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, self.ncol);

        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let nrow = row_offsets.len();
            let ncol = col_offsets.len();
            let mut new_block_polys = vec![vec![T::zero(&self.params); ncol]; nrow];
            let self_block_polys = self.block_entries(row_offsets, col_offsets);
            for i in 0..nrow {
                for j in 0..ncol {
                    new_block_polys[i][j] = self_block_polys[i][j].clone() * rhs;
                }
            }
            new_block_polys
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Neg for MmapMatrix<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let nrow = row_offsets.len();
            let ncol = col_offsets.len();
            let self_block_polys = self.block_entries(row_offsets, col_offsets);
            let mut new_block_polys = vec![vec![T::zero(&self.params); ncol]; nrow];
            for i in 0..nrow {
                for j in 0..ncol {
                    new_block_polys[i][j] = -self_block_polys[i][j].clone();
                }
            }
            new_block_polys
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

fn map_file(file: &File, offset: usize, len: usize) -> Mmap {
    unsafe {
        MmapOptions::new()
            .offset(offset as u64)
            .len(len)
            .populate()
            .map(file)
            .expect("failed to map file")
    }
}

unsafe fn map_file_mut(file: &File, offset: usize, len: usize) -> MmapMut {
    unsafe {
        MmapOptions::new()
            .offset(offset as u64)
            .len(len)
            .populate()
            .map_mut(file)
            .expect("failed to map file")
    }
}

pub fn block_size() -> usize {
    env::var("BLOCK_SIZE").map(|str| str.parse::<usize>().unwrap()).unwrap_or(1000)
}

pub fn block_offsets(rows: Range<usize>, cols: Range<usize>) -> (Vec<usize>, Vec<usize>) {
    let block_size = block_size();
    // *BLOCK_SIZE.get().unwrap();
    let nrow = rows.end - rows.start;
    let ncol = cols.end - cols.start;
    let num_blocks_row = nrow.div_ceil(block_size);
    let num_blocks_col = ncol.div_ceil(block_size);
    let mut row_offsets = vec![rows.start];
    for _ in 0..num_blocks_row {
        let last_row_offset = row_offsets.last().unwrap();
        let sub = rows.end - last_row_offset;
        let len = if sub < block_size { sub } else { block_size };
        row_offsets.push(last_row_offset + len);
    }
    let mut col_offsets = vec![cols.start];
    for _ in 0..num_blocks_col {
        let last_col_offset = col_offsets.last().unwrap();
        let sub = cols.end - last_col_offset;
        let len = if sub < block_size { sub } else { block_size };
        col_offsets.push(last_col_offset + len);
    }
    (row_offsets, col_offsets)
}

fn add_block_matrices<T: MmapMatrixElem>(lhs: Vec<Vec<T>>, rhs: &[Vec<T>]) -> Vec<Vec<T>> {
    let nrow = lhs.len();
    let ncol = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol).map(|j| lhs[i][j].clone() + &rhs[i][j]).collect::<Vec<T>>()
        })
        .collect::<Vec<Vec<T>>>()
}

fn sub_block_matrices<T: MmapMatrixElem>(lhs: Vec<Vec<T>>, rhs: &[Vec<T>]) -> Vec<Vec<T>> {
    let nrow = lhs.len();
    let ncol = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol).map(|j| lhs[i][j].clone() - &rhs[i][j]).collect::<Vec<T>>()
        })
        .collect::<Vec<Vec<T>>>()
}

fn mul_block_matrices<T: MmapMatrixElem>(
    params: &T::Params,
    lhs: Vec<Vec<T>>,
    rhs: Vec<Vec<T>>,
) -> Vec<Vec<T>> {
    let nrow = lhs.len();
    let ncol = rhs[0].len();
    let n_inner = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol)
                .map(|j: usize| {
                    let mut sum = T::zero(params);
                    for k in 0..n_inner {
                        sum += lhs[i][k].clone() * &rhs[k][j];
                    }
                    sum
                })
                .collect::<Vec<T>>()
        })
        .collect::<Vec<Vec<T>>>()
}

// #[cfg(test)]
// #[cfg(feature = "test")]
// mod tests {
//     use std::sync::Arc;

//     use num_bigint::BigUint;

//     use super::*;
//     use crate::poly::{
//         dcrt::{DCRTPolyParams, DCRTPolyUniformSampler},
//         sampler::PolyUniformSampler,
//     };

//     #[test]
//     fn test_matrix_gadget_matrix() {
//         let params = DCRTPolyParams::default();
//         let size = 3;
//         let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, size);
//         assert_eq!(gadget_matrix.size().0, size);
//         assert_eq!(gadget_matrix.size().1, size * params.modulus_bits());
//     }

//     #[test]
//     fn test_matrix_decompose() {
//         let params = DCRTPolyParams::default();
//         let bit_length = params.modulus_bits();

//         // Create a simple 2x8 matrix with some non-zero values
//         let mut matrix_vec = Vec::with_capacity(2);
//         let value = FinRingElem::new(5u32, params.modulus());

//         // Create first row
//         let mut row1 = Vec::with_capacity(8);
//         row1.push(DCRTPoly::from_const(&params, &value));
//         for _ in 1..8 {
//             row1.push(DCRTPoly::const_zero(&params));
//         }

//         // Create second row
//         let mut row2 = Vec::with_capacity(8);
//         row2.push(DCRTPoly::const_zero(&params));
//         row2.push(DCRTPoly::from_const(&params, &value));
//         for _ in 2..8 {
//             row2.push(DCRTPoly::const_zero(&params));
//         }

//         matrix_vec.push(row1);
//         matrix_vec.push(row2);

//         let matrix = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
//         assert_eq!(matrix.size().0, 2);
//         assert_eq!(matrix.size().1, 8);

//         let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, 2);
//         assert_eq!(gadget_matrix.size().0, 2);
//         assert_eq!(gadget_matrix.size().1, 2 * bit_length);

//         let decomposed = matrix.decompose();
//         assert_eq!(decomposed.size().0, 2 * bit_length);
//         assert_eq!(decomposed.size().1, 8);

//         let expected_matrix = gadget_matrix * decomposed;
//         assert_eq!(expected_matrix.size().0, 2);
//         assert_eq!(expected_matrix.size().1, 8);
//         assert_eq!(matrix, expected_matrix);
//     }

//     #[test]
//     fn test_matrix_basic_operations() {
//         let params = DCRTPolyParams::default();

//         // Test zero and identity matrices
//         let zero = DCRTPolyMatrix::zero(&params, 2, 2);
//         let identity = DCRTPolyMatrix::identity(&params, 2, None);

//         // Test matrix creation and equality
//         let value = FinRingElem::new(5u32, params.modulus());

//         // Create a 2x2 matrix with values at (0,0) and (1,1)
//         let matrix_vec = vec![
//             vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
//             vec![DCRTPoly::const_zero(&params), DCRTPoly::from_const(&params, &value)],
//         ];

//         let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
//         assert_eq!(matrix1.entry(0, 0).coeffs()[0], value);
//         let matrix2 = matrix1.clone();
//         assert_eq!(matrix1, matrix2);

//         // Test addition
//         let sum = matrix1.clone() + &matrix2;
//         let value_10 = FinRingElem::new(10u32, params.modulus());
//         assert_eq!(sum.entry(0, 0).coeffs()[0], value_10);

//         // Test subtraction
//         let diff = matrix1.clone() - &matrix2;
//         assert_eq!(diff, zero);

//         // Test multiplication
//         let prod = matrix1 * &identity;
//         assert_eq!(prod.size(), (2, 2));
//         // Check that the product has the same values as the original matrix
//         assert_eq!(prod.entry(0, 0).coeffs()[0], value);
//         assert_eq!(prod.entry(1, 1).coeffs()[0], value);
//     }

//     #[test]
//     fn test_matrix_concatenation() {
//         let params = DCRTPolyParams::default();
//         let value = FinRingElem::new(5u32, params.modulus());

//         // Create first matrix with value at (0,0)
//         let matrix1_vec = vec![
//             vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
//             vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
//         ];

//         let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix1_vec);

//         // Create second matrix with value at (1,1)
//         let matrix2_vec = vec![
//             vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
//             vec![DCRTPoly::const_zero(&params), DCRTPoly::from_const(&params, &value)],
//         ];

//         let matrix2 = DCRTPolyMatrix::from_poly_vec(&params, matrix2_vec);

//         // Test column concatenation
//         let col_concat = matrix1.concat_columns(&[&matrix2]);
//         assert_eq!(col_concat.size().0, 2);
//         assert_eq!(col_concat.size().1, 4);
//         assert_eq!(col_concat.entry(0, 0).coeffs()[0], value);
//         assert_eq!(col_concat.entry(1, 3).coeffs()[0], value);

//         // Test row concatenation
//         let row_concat = matrix1.concat_rows(&[&matrix2]);
//         assert_eq!(row_concat.size().0, 4);
//         assert_eq!(row_concat.size().1, 2);
//         assert_eq!(row_concat.entry(0, 0).coeffs()[0], value);
//         assert_eq!(row_concat.entry(3, 1).coeffs()[0], value);

//         // Test diagonal concatenation
//         let diag_concat = matrix1.concat_diag(&[&matrix2]);
//         assert_eq!(diag_concat.size().0, 4);
//         assert_eq!(diag_concat.size().1, 4);
//         assert_eq!(diag_concat.entry(0, 0).coeffs()[0], value);
//         assert_eq!(diag_concat.entry(3, 3).coeffs()[0], value);
//     }

//     #[test]
//     fn test_matrix_tensor_product() {
//         let params = DCRTPolyParams::default();
//         let value = FinRingElem::new(5u32, params.modulus());

//         // Create first matrix with value at (0,0)
//         let matrix1_vec = vec![
//             vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
//             vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
//         ];

//         let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix1_vec);

//         // Create second matrix with value at (0,0)
//         let matrix2_vec = vec![
//             vec![DCRTPoly::from_const(&params, &value), DCRTPoly::const_zero(&params)],
//             vec![DCRTPoly::const_zero(&params), DCRTPoly::const_zero(&params)],
//         ];

//         let matrix2 = DCRTPolyMatrix::from_poly_vec(&params, matrix2_vec);

//         let tensor = matrix1.tensor(&matrix2);
//         assert_eq!(tensor.size().0, 4);
//         assert_eq!(tensor.size().1, 4);

//         // Check that the (0,0) element is the product of the (0,0) elements
//         let value_25 = FinRingElem::new(25u32, params.modulus());
//         assert_eq!(tensor.entry(0, 0).coeffs()[0], value_25);
//     }

//     #[test]
//     fn test_matrix_modulus_switch() {
//         let params = DCRTPolyParams::default();

//         let value00 = FinRingElem::new(1023782870921908217643761278891282178u128,
// params.modulus());         let value01 =
// FinRingElem::new(8179012198875468938912873783289218738u128, params.modulus());         let
// value10 = FinRingElem::new(2034903202902173762872163465127672178u128, params.modulus());
//         let value11 = FinRingElem::new(1990091289902891278121564387120912660u128,
// params.modulus());

//         let matrix_vec = vec![
//             vec![DCRTPoly::from_const(&params, &value00), DCRTPoly::from_const(&params,
// &value01)],             vec![DCRTPoly::from_const(&params, &value10),
// DCRTPoly::from_const(&params, &value11)],         ];

//         let matrix = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
//         let new_modulus = Arc::new(BigUint::from(2u32));
//         let switched = matrix.modulus_switch(&new_modulus);

//         // Although the value becomes less than the new modulus, the set modulus is still the
// same         assert_eq!(switched.params.modulus(), params.modulus());

//         let new_value00 = value00.modulus_switch(new_modulus.clone());
//         let new_value01 = value01.modulus_switch(new_modulus.clone());
//         let new_value10 = value10.modulus_switch(new_modulus.clone());
//         let new_value11 = value11.modulus_switch(new_modulus.clone());

//         let expected_vec = vec![
//             vec![
//                 DCRTPoly::from_const(&params, &new_value00),
//                 DCRTPoly::from_const(&params, &new_value01),
//             ],
//             vec![
//                 DCRTPoly::from_const(&params, &new_value10),
//                 DCRTPoly::from_const(&params, &new_value11),
//             ],
//         ];

//         let expected = DCRTPolyMatrix::from_poly_vec(&params, expected_vec);
//         assert_eq!(switched, expected);
//     }

//     #[test]
//     #[should_panic(expected = "Addition requires matrices of same dimensions")]
//     #[cfg(debug_assertions)]
//     fn test_matrix_addition_mismatch() {
//         let params = DCRTPolyParams::default();
//         let matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
//         let matrix2 = DCRTPolyMatrix::zero(&params, 2, 3);
//         let _sum = matrix1 + matrix2;
//     }

//     #[test]
//     #[should_panic(expected = "Multiplication condition failed")]
//     #[cfg(debug_assertions)]
//     fn test_matrix_multiplication_mismatch() {
//         let params = DCRTPolyParams::default();
//         let matrix1 = DCRTPolyMatrix::zero(&params, 2, 2);
//         let matrix2 = DCRTPolyMatrix::zero(&params, 3, 2);
//         let _prod = matrix1 * matrix2;
//     }

//     #[test]
//     fn test_matrix_mul_tensor_identity_simple() {
//         let params = DCRTPolyParams::default();
//         let sampler = DCRTPolyUniformSampler::new();

//         // Create matrix S (2x20)
//         let s = sampler.sample_uniform(&params, 2, 20,
// crate::poly::sampler::DistType::FinRingDist);         // Create 'other' matrix (5x7)
//         let other =
//             sampler.sample_uniform(&params, 5, 7, crate::poly::sampler::DistType::FinRingDist);
//         // Perform S * (I_4 ⊗ other)
//         let result = s.mul_tensor_identity(&other, 4);

//         // Check dimensions
//         assert_eq!(result.size().0, 2);
//         assert_eq!(result.size().1, 28);

//         let identity = DCRTPolyMatrix::identity(&params, 4, None);
//         // Check result
//         let expected_result = s * (identity.tensor(&other));

//         assert_eq!(expected_result.size().0, 2);
//         assert_eq!(expected_result.size().1, 28);
//         assert_eq!(result, expected_result)
//     }

//     #[test]
//     fn test_matrix_mul_tensor_identity_decompose_naive() {
//         let params = DCRTPolyParams::default();
//         let sampler = DCRTPolyUniformSampler::new();

//         // Create matrix S (2x2516)
//         let s =
//             sampler.sample_uniform(&params, 2, 2516,
// crate::poly::sampler::DistType::FinRingDist);

//         // Create 'other' matrix (2x13)
//         let other =
//             sampler.sample_uniform(&params, 2, 13, crate::poly::sampler::DistType::FinRingDist);

//         // Decompose 'other' matrix
//         let other_decompose = other.decompose();
//         // Perform S * (I_37 ⊗ G^-1(other))
//         let result: DCRTPolyMatrix = s.mul_tensor_identity(&other_decompose, 37);
//         // Check dimensions
//         assert_eq!(result.size().0, 2);
//         assert_eq!(result.size().1, 481);

//         // Check result
//         let tensor = identity_tensor_matrix(37, &other_decompose);
//         let expected_result = s * tensor;

//         assert_eq!(expected_result.size().0, 2);
//         assert_eq!(expected_result.size().1, 481);
//         assert_eq!(result, expected_result)
//     }

//     #[test]
//     fn test_matrix_mul_tensor_identity_decompose_optimal() {
//         let params = DCRTPolyParams::default();
//         let sampler = DCRTPolyUniformSampler::new();

//         // Create matrix S (2x2516)
//         let s =
//             sampler.sample_uniform(&params, 2, 2516,
// crate::poly::sampler::DistType::FinRingDist);

//         // Create 'other' matrix (2x13)
//         let other =
//             sampler.sample_uniform(&params, 2, 13, crate::poly::sampler::DistType::FinRingDist);

//         // Perform S * (I_37 ⊗ G^-1(other))
//         let result: DCRTPolyMatrix = s.mul_tensor_identity_decompose(&other, 37);

//         // Check dimensions
//         assert_eq!(result.size().0, 2);
//         assert_eq!(result.size().1, 481);

//         // Check result
//         let decomposed = other.decompose();
//         let tensor = identity_tensor_matrix(37, &decomposed);
//         let expected_result_1 = s.clone() * tensor;
//         let expected_result_2 = s.mul_tensor_identity(&decomposed, 37);
//         assert_eq!(expected_result_1, expected_result_2);

//         assert_eq!(expected_result_1.size().0, 2);
//         assert_eq!(expected_result_1.size().1, 481);

//         assert_eq!(expected_result_2.size().0, 2);
//         assert_eq!(expected_result_2.size().1, 481);

//         assert_eq!(result, expected_result_1);
//         assert_eq!(result, expected_result_2);
//     }

//     fn identity_tensor_matrix(identity_size: usize, matrix: &DCRTPolyMatrix) -> DCRTPolyMatrix {
//         let mut others = vec![];
//         for _ in 1..identity_size {
//             others.push(matrix);
//         }
//         matrix.concat_diag(&others[..])
//     }
// }
