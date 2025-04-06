use crate::{parallel_iter, utils::debug_mem};
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

    pub fn block_entries_column(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        let entry_size = self.entry_size();
        let total_cols = self.ncol;
        let nrows = rows.end - rows.start;
        let start_offset = entry_size * (rows.start * total_cols);
        let total_bytes = entry_size * (nrows * total_cols);
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) as usize };
        let aligned_offset = start_offset - (start_offset % page_size);
        let offset_in_page = start_offset - aligned_offset;
        let mapping_len = ((offset_in_page + total_bytes + page_size - 1) / page_size) * page_size;
        let mmap = map_file(&self.file, aligned_offset, mapping_len);
        let data = mmap[offset_in_page..offset_in_page + total_bytes].to_vec();
        drop(mmap);
        let columns: Vec<Vec<T>> = (cols.start..cols.end)
            .into_par_iter()
            .map(|col| {
                let mut col_vec = Vec::with_capacity(nrows);
                for row in 0..nrows {
                    let offset = row * total_cols * entry_size + col * entry_size;
                    let elem_bytes = &data[offset..offset + entry_size];
                    col_vec.push(T::from_bytes_to_elem(&self.params, elem_bytes));
                }
                col_vec
            })
            .collect();

        columns
    }

    pub fn block_entries_row(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        let entry_size = self.entry_size();
        let num_cols = cols.end - cols.start;
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) as usize };

        parallel_iter!(rows)
            .map(|i| {
                let desired_offset = entry_size * (i * self.ncol + cols.start);
                let desired_len = entry_size * num_cols;
                let aligned_offset = desired_offset - (desired_offset % page_size);
                let offset_in_page = desired_offset - aligned_offset;
                let mapping_len =
                    ((offset_in_page + desired_len + page_size - 1) / page_size) * page_size;
                let mmap = map_file(&self.file, aligned_offset, mapping_len);
                let row_bytes = mmap[offset_in_page..offset_in_page + desired_len].to_vec();
                let row_entries = row_bytes
                    .chunks(entry_size)
                    .map(|entry| T::from_bytes_to_elem(&self.params, entry))
                    .collect_vec();
                drop(mmap);
                row_entries
            })
            .collect()
    }

    /// # Output Interface Distinction
    /// **Important:** This function differs from `replace_entries_row` in the expected layout of
    /// the output from the closure `f`. For `replace_entries_column`, the returned 2D vector
    /// should be arranged in a column-major format:
    /// - The **first dimension** (outer vector) corresponds to the block’s **columns**.
    /// - The **second dimension** (inner vector) corresponds to the block’s **rows**.
    /// In other words, the first and second dimensions correspond to the column and row,
    /// respectively.
    pub fn replace_entries_column<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let (row_offsets, col_offsets) = block_offsets(rows, cols);
        parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_col_idx, next_block_col_idx)| {
                parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_row_idx, next_block_row_idx)| {
                        let new_entries = f(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // The non-overlapping nature of blocks ensures thread safety.
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

    pub fn replace_entries_row<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let (row_offsets, col_offsets) = block_offsets(rows, cols);
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
        let (offsets, _) = block_offsets(diags, 0..0);
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
        let block_size = block_size();
        let row_block_size = block_size.div_ceil(row_scale);
        let col_block_size = block_size.div_ceil(col_scale);
        let (row_offsets, col_offsets) =
            block_offsets_distinct_block_sizes(rows, cols, row_block_size, col_block_size);
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
                                *cur_block_row_idx * row_scale..*next_block_row_idx * row_scale,
                                *cur_block_col_idx * col_scale..*next_block_col_idx * col_scale,
                                new_entries,
                            );
                        }
                    },
                );
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
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) as usize };
        parallel_iter!(rows).for_each(|i| {
            let desired_offset = entry_size * (i * self.ncol + cols.start);
            let aligned_offset = desired_offset - (desired_offset % page_size);
            let offset_in_page = desired_offset - aligned_offset;
            let desired_len = entry_size * cols.len();
            let mapping_len = {
                // Calculate the total bytes needed, then round up to a multiple of page_size.
                let total = offset_in_page + desired_len;
                (total + page_size - 1) / page_size * page_size
            };
            let mut mmap = unsafe { map_file_mut(&self.file, aligned_offset, mapping_len) };
            let region = &mut mmap[offset_in_page..offset_in_page + desired_len];
            let bytes = new_entries[i - row_start]
                .iter()
                .flat_map(|poly| poly.as_elem_to_bytes())
                .collect::<Vec<_>>();
            region.copy_from_slice(&bytes);
            drop(mmap);
        });
    }

    pub fn entry(&self, i: usize, j: usize) -> T {
        self.block_entries_row(i..i + 1, j..j + 1)[0][0].clone()
    }

    pub fn get_row(&self, i: usize) -> Vec<T> {
        self.block_entries_row(i..i + 1, 0..self.ncol)[0].clone()
    }

    pub fn get_column(&self, j: usize) -> Vec<T> {
        self.block_entries_column(0..self.nrow, j..j + 1)[0].clone()
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
            self.block_entries_row(row_offsets, col_offsets)
        };
        new_matrix.replace_entries_row(0..nrow, 0..ncol, f);
        new_matrix
    }

    pub fn zero(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        Self::new_empty(params, nrow, ncol)
    }

    pub fn identity(params: &T::Params, size: usize, scalar: Option<T>) -> Self {
        let mut new_matrix = Self::new_empty(params, size, size);
        let scalar = scalar.unwrap_or_else(|| T::one(params));
        let f = |offsets: Range<usize>| -> Vec<Vec<T>> {
            let len = offsets.len();
            (0..len)
                .map(|i| {
                    (0..len)
                        .map(|j| if i == j { scalar.clone() } else { T::zero(params) })
                        .collect()
                })
                .collect()
        };
        new_matrix.replace_entries_diag(0..size, f);
        new_matrix
    }

    pub fn transpose(&self) -> Self {
        let mut new_matrix = Self::new_empty(&self.params, self.ncol, self.nrow);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let cur_entries = self.block_entries_row(col_offsets.clone(), row_offsets.clone());
            let row_offsets_len = row_offsets.len();
            let col_offsets_len = col_offsets.len();
            (0..row_offsets_len)
                .into_par_iter()
                .map(|i| {
                    (0..col_offsets_len)
                        .into_par_iter()
                        .map(|j| cur_entries[j][i].clone())
                        .collect::<Vec<T>>()
                })
                .collect::<Vec<Vec<T>>>()
        };
        new_matrix.replace_entries_row(0..self.ncol, 0..self.nrow, f);
        new_matrix
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    pub fn concat_columns(&self, others: &[&Self]) -> Self {
        #[cfg(debug_assertions)]
        for (idx, other) in others.iter().enumerate() {
            if self.nrow != other.nrow {
                panic!(
                    "Concat error: while the shape of the first matrix is ({}, {}), \
                     that of the {}-th matrix is ({},{})",
                    self.nrow, self.ncol, idx, other.nrow, other.ncol
                );
            }
        }
        let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);
        let mut combined_data: Vec<Vec<T>> = Vec::with_capacity(updated_ncol);
        combined_data.extend(self.block_entries_column(0..self.nrow, 0..self.ncol));
        for other in others {
            combined_data.extend(other.block_entries_column(0..self.nrow, 0..other.ncol));
        }
        let mut new_matrix = Self::new_empty(&self.params, self.nrow, updated_ncol);
        new_matrix.replace_entries_column(0..self.nrow, 0..updated_ncol, |row_range, col_range| {
            row_range
                .map(|r| {
                    (col_range.start..col_range.end)
                        .map(|c| combined_data[c][r].clone())
                        .collect::<Vec<T>>()
                })
                .collect::<Vec<Vec<T>>>()
        });
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
        let mut combined_data: Vec<Vec<T>> = Vec::with_capacity(updated_nrow);
        combined_data.extend(self.block_entries_row(0..self.nrow, 0..self.ncol));
        for other in others {
            combined_data.extend(other.block_entries_row(0..other.nrow, 0..self.ncol));
        }
        let mut new_matrix = Self::new_empty(&self.params, updated_nrow, self.ncol);
        new_matrix.replace_entries_row(0..updated_nrow, 0..self.ncol, |row_range, col_range| {
            row_range.map(|r| combined_data[r][col_range.clone()].to_vec()).collect()
        });

        new_matrix
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    pub fn concat_diag(&self, others: &[&Self]) -> Self {
        let updated_nrow = others.iter().fold(self.nrow, |acc, other| acc + other.nrow);
        let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);

        let mut new_matrix = Self::new_empty(&self.params, updated_nrow, updated_ncol);
        let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries_row(row_offsets, col_offsets)
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..self.ncol, self_f);
        debug_mem("self replaced in concat_diag");

        let mut row_acc = self.nrow;
        let mut col_acc = self.ncol;
        for (idx, other) in others.iter().enumerate() {
            let other_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                let row_offsets = row_offsets.start - row_acc..row_offsets.end - row_acc;
                let col_offsets = col_offsets.start - col_acc..col_offsets.end - col_acc;
                other.block_entries_row(row_offsets, col_offsets)
            };
            new_matrix.replace_entries_row(
                row_acc..row_acc + other.nrow,
                col_acc..col_acc + other.ncol,
                other_f,
            );
            debug_mem(format!("the {}-th other replaced in concat_diag", idx));
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
                let sub_matrix = other * scalar;
                parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_row_idx, next_block_row_idx)| {
                        parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                            |(cur_block_col_idx, next_block_col_idx)| {
                                let sub_block_polys = sub_matrix.block_entries_row(
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
            self.block_entries_row(row_offsets, col_offsets)
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..self.ncol, f);
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
                        let self_block_polys = self.block_entries_row(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let other_block_polys = other.block_entries_row(
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
        // debug_mem("Drop MmapMatrix");
        self.file.set_len(0).expect("failed to truncate file");
        // debug_mem("Truncate file");
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
            let self_block_polys = self.block_entries_row(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries_row(row_offsets, col_offsets);
            add_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Sub<MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<T: MmapMatrixElem> Sub<&MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn sub(self, rhs: &Self) -> Self::Output {
        &self - rhs
    }
}

impl<T: MmapMatrixElem> Sub<&MmapMatrix<T>> for &MmapMatrix<T> {
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
            let self_block_polys = self.block_entries_row(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries_row(row_offsets, col_offsets);
            sub_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Mul<MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: MmapMatrixElem> Mul<&MmapMatrix<T>> for MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<T: MmapMatrixElem> Mul<&MmapMatrix<T>> for &MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: &MmapMatrix<T>) -> Self::Output {
        debug_assert!(
            self.ncol == rhs.nrow,
            "Multiplication condition failed: self.ncol ({}) must equal rhs.nrow ({})",
            self.ncol,
            rhs.nrow
        );
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, rhs.ncol);
        let (_, ip_offsets) = block_offsets(0..0, 0..self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            ip_offsets
                .iter()
                .tuple_windows()
                .map(|(cur_block_ip_idx, next_block_ip_idx)| {
                    let self_block_polys = self.block_entries_row(
                        row_offsets.clone(),
                        *cur_block_ip_idx..*next_block_ip_idx,
                    );
                    let other_block_polys = rhs.block_entries_row(
                        *cur_block_ip_idx..*next_block_ip_idx,
                        col_offsets.clone(),
                    );
                    mul_block_matrices(self_block_polys, other_block_polys)
                })
                .reduce(|acc, muled| add_block_matrices(muled, &acc))
                .unwrap()
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..rhs.ncol, f);
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
        &self * rhs
    }
}

impl<T: MmapMatrixElem> Mul<T> for &MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MmapMatrixElem> Mul<&T> for &MmapMatrix<T> {
    type Output = MmapMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries_row(row_offsets, col_offsets)
                .into_iter()
                .map(|row| row.into_iter().map(|elem| elem * rhs).collect::<Vec<T>>())
                .collect::<Vec<Vec<T>>>()
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MmapMatrixElem> Neg for MmapMatrix<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut new_matrix = MmapMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries_row(row_offsets, col_offsets)
                .into_iter()
                .map(|row| row.into_iter().map(|elem| -elem.clone()).collect())
                .collect()
        };
        new_matrix.replace_entries_row(0..self.nrow, 0..self.ncol, f);
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
    env::var("BLOCK_SIZE").map(|str| str.parse::<usize>().unwrap()).unwrap_or(100)
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
    // *BLOCK_SIZE.get().unwrap();
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

fn mul_block_matrices<T: MmapMatrixElem>(lhs: Vec<Vec<T>>, rhs: Vec<Vec<T>>) -> Vec<Vec<T>> {
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
