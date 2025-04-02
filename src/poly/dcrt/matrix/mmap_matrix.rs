use crate::parallel_iter;
use itertools::Itertools;
use memmap2::{Mmap, MmapMut, MmapOptions};
// use once_cell::sync::OnceCell;
use rayon::prelude::*;
use std::{
    cell::RefCell,
    env,
    fmt::Debug,
    fs::File,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Range, Sub, SubAssign},
    sync::{mpsc::channel, Arc, RwLock},
    thread,
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

// Thread-local buffer to store working data
struct ThreadLocalBuffer<T: MmapMatrixElem> {
    buffer: Arc<RwLock<Vec<Vec<T>>>>,
}

impl<T: MmapMatrixElem> ThreadLocalBuffer<T> {
    fn new() -> Self {
        Self { buffer: Arc::new(RwLock::new(Vec::new())) }
    }

    fn resize(&self, rows: usize, cols: usize, params: &T::Params) {
        let mut buffer = self.buffer.write().unwrap();
        buffer.resize(rows, Vec::with_capacity(cols));
        for row in buffer.iter_mut() {
            row.resize(cols, T::zero(params));
        }
    }

    fn get(&self) -> &Arc<RwLock<Vec<Vec<T>>>> {
        &self.buffer
    }
}

impl<T: MmapMatrixElem> Clone for ThreadLocalBuffer<T> {
    fn clone(&self) -> Self {
        Self { buffer: self.buffer.clone() }
    }
}

pub struct MmapMatrix<T: MmapMatrixElem> {
    pub params: T::Params,
    pub file: File,
    pub nrow: usize,
    pub ncol: usize,
    buffer: ThreadLocalBuffer<T>,
}

impl<T: MmapMatrixElem> MmapMatrix<T> {
    pub fn new_empty(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        let entry_size = params.entry_size();
        let len = entry_size * nrow * ncol;
        let file = tempfile().expect("failed to open file");
        file.set_len(len as u64).expect("failed to set file length");
        Self { params: params.clone(), file, nrow, ncol, buffer: ThreadLocalBuffer::new() }
    }

    pub fn entry_size(&self) -> usize {
        self.params.entry_size()
    }

    pub fn block_entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        let entry_size = self.entry_size();
        let nrows = rows.len();
        let ncols = cols.len();

        // Resize thread-local buffer if needed
        self.buffer.resize(nrows, ncols, &self.params);

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
        let entry_size = self.entry_size();
        let ncol = self.ncol;

        let (tx, rx) = channel::<(usize, Vec<u8>)>();
        let file = self.file.try_clone().expect("failed to clone file");

        let io_handle = thread::spawn(move || {
            while let Ok((offset, data)) = rx.recv() {
                unsafe {
                    let mut mmap = map_file_mut(&file, offset, data.len());
                    mmap.copy_from_slice(&data);
                    drop(mmap);
                }
            }
        });

        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        // Compute new entries in thread-local buffer
                        let new_entries = f(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );

                        // Write back to mmap through I/O thread
                        send_block_rows(
                            &tx,
                            new_entries,
                            *cur_block_row_idx,
                            *cur_block_col_idx,
                            ncol,
                            entry_size,
                        );
                    },
                );
            },
        );
        drop(tx);
        io_handle.join().expect("I/O thread panicked");
    }

    /// Replace diagonal entries using batching via a dedicated I/O thread.
    pub fn replace_entries_diag<F>(&mut self, diags: Range<usize>, f: F)
    where
        F: Fn(Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        let (offsets, _) = block_offsets(diags.clone(), 0..0);
        let entry_size = self.entry_size();
        let ncol = self.ncol;

        let (tx, rx) = channel::<(usize, Vec<u8>)>();
        let file = self.file.try_clone().expect("failed to clone file");

        let io_handle = thread::spawn(move || {
            while let Ok((offset, data)) = rx.recv() {
                unsafe {
                    let mut mmap = map_file_mut(&file, offset, data.len());
                    mmap.copy_from_slice(&data);
                    drop(mmap);
                }
            }
        });

        parallel_iter!(offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_diag_idx, next_block_diag_idx)| {
                let new_entries = f(*cur_block_diag_idx..*next_block_diag_idx);
                // Here, rows and cols are the same for diagonal blocks.
                send_block_rows(
                    &tx,
                    new_entries,
                    *cur_block_diag_idx,
                    *cur_block_diag_idx,
                    ncol,
                    entry_size,
                );
            },
        );

        drop(tx);
        io_handle.join().expect("I/O thread panicked");
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
            let cur_entries = self.block_entries(col_offsets.clone(), row_offsets.clone());
            row_offsets
                .map(|i| col_offsets.clone().map(|j| cur_entries[j][i].clone()).collect::<Vec<T>>())
                .collect::<Vec<Vec<T>>>()
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
        let (tx, rx) = channel::<(usize, Vec<u8>)>();
        let file = new_matrix.file.try_clone().expect("failed to clone file");
        let entry_size = new_matrix.entry_size();
        let ncol = new_matrix.ncol;
        let io_handle = thread::spawn(move || {
            while let Ok((offset, data)) = rx.recv() {
                unsafe {
                    let mut mmap = map_file_mut(&file, offset, data.len());
                    mmap.copy_from_slice(&data);
                    drop(mmap);
                }
            }
        });

        parallel_iter!(0..self.nrow).for_each(|i| {
            parallel_iter!(0..self.ncol).for_each(|j| {
                let scalar = self.entry(i, j);
                let sub_matrix = other * scalar;
                parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_row_idx, next_block_row_idx)| {
                        parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                            |(cur_block_col_idx, next_block_col_idx)| {
                                let sub_block_polys = sub_matrix.block_entries(
                                    *cur_block_row_idx..*next_block_row_idx,
                                    *cur_block_col_idx..*next_block_col_idx,
                                );
                                let global_row = i * sub_matrix.nrow + *cur_block_row_idx;
                                let global_col = j * sub_matrix.ncol + *cur_block_col_idx;
                                send_block_rows(
                                    &tx,
                                    sub_block_polys,
                                    global_row,
                                    global_col,
                                    ncol,
                                    entry_size,
                                );
                            },
                        );
                    },
                );
            });
        });
        drop(tx);
        io_handle.join().expect("I/O thread panicked");
        new_matrix
    }

    fn add_block_matrices(&self, lhs: Vec<Vec<T>>, rhs: &[Vec<T>]) -> Vec<Vec<T>> {
        let nrow = lhs.len();
        let ncol = lhs[0].len();

        self.buffer.resize(nrow, ncol, &self.params);
        let buffer = self.buffer.get().clone();

        let result: Vec<Vec<T>> = parallel_iter!(0..nrow)
            .map(|i| parallel_iter!(0..ncol).map(|j| lhs[i][j].clone() + &rhs[i][j]).collect())
            .collect();

        *buffer.write().unwrap() = result.clone();
        result
    }

    fn sub_block_matrices(&self, lhs: Vec<Vec<T>>, rhs: &[Vec<T>]) -> Vec<Vec<T>> {
        let nrow = lhs.len();
        let ncol = lhs[0].len();

        self.buffer.resize(nrow, ncol, &self.params);
        let buffer = self.buffer.get().clone();

        let result: Vec<Vec<T>> = parallel_iter!(0..nrow)
            .map(|i| parallel_iter!(0..ncol).map(|j| lhs[i][j].clone() - &rhs[i][j]).collect())
            .collect();

        *buffer.write().unwrap() = result.clone();
        result
    }

    fn mul_block_matrices(&self, lhs: Vec<Vec<T>>, rhs: Vec<Vec<T>>) -> Vec<Vec<T>> {
        let nrow = lhs.len();
        let ncol = rhs[0].len();
        let n_inner = lhs[0].len();

        self.buffer.resize(nrow, ncol, &self.params);
        let buffer = self.buffer.get().clone();

        let result: Vec<Vec<T>> = parallel_iter!(0..nrow)
            .map(|i| {
                parallel_iter!(0..ncol)
                    .map(|j| {
                        (0..n_inner)
                            .map(|k| lhs[i][k].clone() * &rhs[k][j])
                            .reduce(|acc, prod| acc + prod)
                            .unwrap()
                    })
                    .collect()
            })
            .collect();

        *buffer.write().unwrap() = result.clone();
        result
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
            self.add_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
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
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            self.sub_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
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
                    let self_block_polys = self
                        .block_entries(row_offsets.clone(), *cur_block_ip_idx..*next_block_ip_idx);
                    let other_block_polys = rhs
                        .block_entries(*cur_block_ip_idx..*next_block_ip_idx, col_offsets.clone());
                    self.mul_block_matrices(self_block_polys, other_block_polys)
                })
                .reduce(|acc, muled| self.add_block_matrices(muled, &acc))
                .unwrap()
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
            self.block_entries(row_offsets, col_offsets)
                .into_iter()
                .map(|row| row.into_iter().map(|elem| elem * rhs).collect::<Vec<T>>())
                .collect::<Vec<Vec<T>>>()
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
            self.block_entries(row_offsets, col_offsets)
                .into_iter()
                .map(|row| row.into_iter().map(|elem| -elem.clone()).collect())
                .collect()
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

pub unsafe fn map_file_mut(file: &File, offset: usize, len: usize) -> MmapMut {
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

pub fn send_block_rows<T: MmapMatrixElem>(
    tx: &std::sync::mpsc::Sender<(usize, Vec<u8>)>,
    new_entries: Vec<Vec<T>>,
    row_start: usize,
    col_start: usize,
    ncol: usize,
    entry_size: usize,
) {
    // For each row in the block, flatten its entries to a Vec<u8> and send (offset, data)
    for (i, row) in new_entries.into_iter().enumerate() {
        let global_row = row_start + i;
        let offset = entry_size * (global_row * ncol + col_start);
        let bytes: Vec<u8> = row.into_iter().flat_map(|elem| elem.as_elem_to_bytes()).collect();
        tx.send((offset, bytes)).expect("failed to send I/O batch");
    }
}
