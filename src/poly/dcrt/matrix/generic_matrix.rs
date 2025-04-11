#[cfg(feature = "disk")]
use crate::utils::debug_mem;
use crate::{
    parallel_iter,
    poly::{MatrixElem, MatrixParams},
    utils::block_size,
};
use itertools::Itertools;
#[cfg(feature = "disk")]
use memmap2::{Mmap, MmapMut, MmapOptions};
use rayon::prelude::*;
#[cfg(feature = "disk")]
use std::fs::File;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Range, Sub},
};
// use sysinfo::System;
#[cfg(feature = "disk")]
use libc;
#[cfg(feature = "disk")]
use tempfile::tempfile;

// static BLOCK_SIZE: OnceCell<usize> = OnceCell::new();

#[cfg_attr(not(feature = "disk"), derive(Clone))]
pub struct GenericMatrix<T: MatrixElem> {
    pub params: T::Params,
    #[cfg(feature = "disk")]
    file: File,
    #[cfg(not(feature = "disk"))]
    inner: Vec<Vec<T>>,
    pub nrow: usize,
    pub ncol: usize,
}

impl<T: MatrixElem> GenericMatrix<T> {
    pub fn new_empty(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        #[cfg(feature = "disk")]
        {
            let entry_size = params.entry_size();
            let len = entry_size * nrow * ncol;
            let file = tempfile().expect("failed to open file");
            file.set_len(len as u64).expect("failed to set file length");
            Self { params: params.clone(), file, nrow, ncol }
        }
        #[cfg(not(feature = "disk"))]
        {
            let inner = vec![vec![T::zero(params); ncol]; nrow];
            Self { params: params.clone(), inner, nrow, ncol }
        }
    }

    pub fn entry_size(&self) -> usize {
        self.params.entry_size()
    }

    pub fn block_entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<T>> {
        #[cfg(feature = "disk")]
        let entry_size = self.entry_size();
        #[cfg(feature = "disk")]
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) as usize };
        parallel_iter!(rows)
            .map(|i| {
                #[cfg(feature = "disk")]
                {
                    let raw_offset = entry_size * (i * self.ncol + cols.start);
                    let aligned_offset = raw_offset - (raw_offset % page_size);
                    let offset_adjustment = raw_offset - aligned_offset;
                    let required_size = offset_adjustment + entry_size * cols.len();
                    let mapping_size = if required_size % page_size == 0 {
                        required_size
                    } else {
                        ((required_size / page_size) + 1) * page_size
                    };
                    let mmap = map_file(&self.file, aligned_offset, mapping_size);
                    let row_data =
                        &mmap[offset_adjustment..offset_adjustment + entry_size * cols.len()];
                    let row_col_vec = row_data
                        .chunks(entry_size)
                        .map(|entry| T::from_bytes_to_elem(&self.params, entry))
                        .collect_vec();
                    drop(mmap);
                    row_col_vec
                }
                #[cfg(not(feature = "disk"))]
                {
                    self.inner[i][cols.start..cols.end].to_vec()
                }
            })
            .collect::<Vec<Vec<_>>>()
    }

    pub fn replace_entries<F>(&mut self, rows: Range<usize>, cols: Range<usize>, f: F)
    where
        F: Fn(Range<usize>, Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        // debug_mem(format!(
        //     "replace_entries: row_offsets: {:?}, col_offsets: {:?}",
        //     row_offsets, col_offsets
        // ));
        if self.nrow == 0 || self.ncol == 0 {
            return;
        }
        #[cfg(feature = "disk")]
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

        #[cfg(not(feature = "disk"))]
        {
            let polys = f(rows.clone(), cols.clone());
            debug_assert_eq!(polys.len(), rows.len());
            debug_assert_eq!(polys[0].len(), cols.len());
            self.inner[rows.start..rows.end].par_iter_mut().enumerate().for_each(
                |(i, row_data)| {
                    row_data[cols.start..cols.end].clone_from_slice(&polys[i]);
                },
            );
            //     let entry_size = self.entry_size();
            //     let mut_mmap = &mut self.mmap
            //         [entry_size * (rows.start * self.ncol)..entry_size * (rows.end * self.ncol)];
            //     mut_mmap.par_chunks_mut(entry_size * self.ncol).enumerate().for_each(
            //         |(i, row_data)| {
            //             let row_mut_mmap =
            //                 &mut row_data[entry_size * (cols.start)..entry_size * (cols.end)];
            //             row_mut_mmap.par_chunks_mut(entry_size).enumerate().for_each(
            //                 |(j, entry_data)| {
            //                     let poly = f(
            //                         rows.start + i..rows.start + i + 1,
            //                         cols.start + j..cols.start + j + 1,
            //                     )[0][0]
            //                         .clone();
            //                     let bytes = poly.as_elem_to_bytes();
            //                     entry_data.copy_from_slice(&bytes);
            //                 },
            //             );
            //         },
            //     );
        }
    }

    pub fn replace_entries_diag<F>(&mut self, diags: Range<usize>, f: F)
    where
        F: Fn(Range<usize>) -> Vec<Vec<T>> + Send + Sync,
    {
        #[cfg(feature = "disk")]
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

        #[cfg(not(feature = "disk"))]
        {
            let polys = f(diags.clone());
            debug_assert_eq!(polys.len(), diags.len());
            debug_assert_eq!(polys[0].len(), diags.len());
            self.inner[diags.start..diags.end].par_iter_mut().enumerate().for_each(
                |(i, row_data)| {
                    row_data[diags.start..diags.end].clone_from_slice(&polys[i]);
                },
            );
        }
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
        #[cfg(feature = "disk")]
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

        #[cfg(not(feature = "disk"))]
        {
            let polys = f(rows.clone(), cols.clone());
            debug_assert_eq!(polys.len(), rows.len() * row_scale);
            debug_assert_eq!(polys[0].len(), cols.len() * col_scale);
            self.inner[rows.start * row_scale..rows.end * row_scale]
                .par_iter_mut()
                .enumerate()
                .for_each(|(i, row_data)| {
                    row_data[cols.start * col_scale..cols.end * col_scale]
                        .clone_from_slice(&polys[i]);
                });
        }
    }

    #[cfg(feature = "disk")]
    pub(crate) unsafe fn replace_block_entries(
        &self,
        rows: Range<usize>,
        cols: Range<usize>,
        new_entries: Vec<Vec<T>>,
    ) {
        debug_assert_eq!(new_entries.len(), rows.end - rows.start);
        debug_assert_eq!(new_entries[0].len(), cols.end - cols.start);
        let entry_size = self.entry_size();
        let num_cols = cols.end - cols.start;
        let row_start = rows.start;
        let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) as usize };
        parallel_iter!(rows).for_each(|i| {
            let desired_offset = entry_size * (i * self.ncol + cols.start);
            let aligned_offset = desired_offset - (desired_offset % page_size);
            let offset_in_page = desired_offset - aligned_offset;
            let required_len = offset_in_page + entry_size * num_cols;
            let mapping_len = required_len.div_ceil(page_size) * page_size;
            let bytes = new_entries[i - row_start]
                .iter()
                .flat_map(|poly| poly.as_elem_to_bytes())
                .collect::<Vec<_>>();
            let mut mmap = unsafe { map_file_mut(&self.file, aligned_offset, mapping_len) };
            mmap[offset_in_page..offset_in_page + entry_size * num_cols].copy_from_slice(&bytes);
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
        #[cfg(feature = "disk")]
        {
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
        #[cfg(not(feature = "disk"))]
        {
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
    }

    pub fn zero(params: &T::Params, nrow: usize, ncol: usize) -> Self {
        Self::new_empty(params, nrow, ncol)
    }

    pub fn identity(params: &T::Params, size: usize, scalar: Option<T>) -> Self {
        #[cfg(feature = "disk")]
        {
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
        #[cfg(not(feature = "disk"))]
        {
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
    }

    pub fn transpose(&self) -> Self {
        #[cfg(feature = "disk")]
        {
            let mut new_matrix = Self::new_empty(&self.params, self.ncol, self.nrow);
            let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                let cur_entries = self.block_entries(col_offsets.clone(), row_offsets.clone());
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
            new_matrix.replace_entries(0..self.ncol, 0..self.nrow, f);
            new_matrix
        }
        #[cfg(not(feature = "disk"))]
        {
            let nrow = self.ncol;
            let ncol = self.nrow;
            let inner: Vec<Vec<T>> = (0..self.ncol)
                .into_par_iter()
                .map(|i| {
                    (0..self.nrow)
                        .into_par_iter()
                        .map(|j| self.inner[j][i].clone())
                        .collect::<Vec<T>>()
                })
                .collect();

            Self { inner, params: self.params.clone(), nrow, ncol }
        }
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
        #[cfg(feature = "disk")]
        {
            let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);
            let mut new_matrix = Self::new_empty(&self.params, self.nrow, updated_ncol);
            let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                self.block_entries(row_offsets, col_offsets)
            };
            new_matrix.replace_entries(0..self.nrow, 0..self.ncol, self_f);
            debug_mem("self replaced in concat_columns");

            let mut col_acc = self.ncol;
            for (idx, other) in others.iter().enumerate() {
                let other_f =
                    |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                        let col_offsets = col_offsets.start - col_acc..col_offsets.end - col_acc;
                        other.block_entries(row_offsets, col_offsets)
                    };
                new_matrix.replace_entries(0..self.nrow, col_acc..col_acc + other.ncol, other_f);
                debug_mem(format!("the {}-th other replaced in concat_columns", idx));
                col_acc += other.ncol;
            }
            debug_assert_eq!(col_acc, updated_ncol);
            new_matrix
        }
        #[cfg(not(feature = "disk"))]
        {
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
        #[cfg(feature = "disk")]
        {
            let updated_nrow = others.iter().fold(self.nrow, |acc, other| acc + other.nrow);

            let mut new_matrix = Self::new_empty(&self.params, updated_nrow, self.ncol);
            let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                self.block_entries(row_offsets, col_offsets)
            };
            new_matrix.replace_entries(0..self.nrow, 0..self.ncol, self_f);
            debug_mem("self replaced in concat_rows");

            let mut row_acc = self.nrow;
            for (idx, other) in others.iter().enumerate() {
                let other_f =
                    |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                        let row_offsets = row_offsets.start - row_acc..row_offsets.end - row_acc;
                        other.block_entries(row_offsets, col_offsets)
                    };
                new_matrix.replace_entries(row_acc..row_acc + other.nrow, 0..self.ncol, other_f);
                debug_mem(format!("the {}-th other replaced in concat_rows", idx));
                row_acc += other.nrow;
            }
            debug_assert_eq!(row_acc, updated_nrow);
            new_matrix
        }
        #[cfg(not(feature = "disk"))]
        {
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
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    pub fn concat_diag(&self, others: &[&Self]) -> Self {
        #[cfg(feature = "disk")]
        {
            let updated_nrow = others.iter().fold(self.nrow, |acc, other| acc + other.nrow);
            let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);

            let mut new_matrix = Self::new_empty(&self.params, updated_nrow, updated_ncol);
            let self_f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                self.block_entries(row_offsets, col_offsets)
            };
            new_matrix.replace_entries(0..self.nrow, 0..self.ncol, self_f);
            debug_mem("self replaced in concat_diag");

            let mut row_acc = self.nrow;
            let mut col_acc = self.ncol;
            for (idx, other) in others.iter().enumerate() {
                let other_f =
                    |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
                        let row_offsets = row_offsets.start - row_acc..row_offsets.end - row_acc;
                        let col_offsets = col_offsets.start - col_acc..col_offsets.end - col_acc;
                        other.block_entries(row_offsets, col_offsets)
                    };
                new_matrix.replace_entries(
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
        #[cfg(not(feature = "disk"))]
        {
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
    }

    pub fn tensor(&self, other: &Self) -> Self {
        #[cfg(feature = "disk")]
        {
            let new_matrix =
                Self::new_empty(&self.params, self.nrow * other.nrow, self.ncol * other.ncol);
            let (row_offsets, col_offsets) = block_offsets(0..other.nrow, 0..other.ncol);
            parallel_iter!(0..self.nrow).for_each(|i| {
                parallel_iter!(0..self.ncol).for_each(|j| {
                    let scalar = self.entry(i, j);
                    let sub_matrix = other * scalar;
                    parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
                        |(cur_block_row_idx, next_block_row_idx)| {
                            parallel_iter!(col_offsets.iter().tuple_windows().collect_vec())
                                .for_each(|(cur_block_col_idx, next_block_col_idx)| {
                                    let sub_block_polys = sub_matrix.block_entries(
                                        *cur_block_row_idx..*next_block_row_idx,
                                        *cur_block_col_idx..*next_block_col_idx,
                                    );
                                    // This is secure because the modified entries are not
                                    // overlapped among threads
                                    unsafe {
                                        new_matrix.replace_block_entries(
                                            i * sub_matrix.nrow + *cur_block_row_idx..
                                                i * sub_matrix.nrow + *next_block_row_idx,
                                            j * sub_matrix.ncol + *cur_block_col_idx..
                                                j * sub_matrix.ncol + *next_block_col_idx,
                                            sub_block_polys,
                                        );
                                    }
                                });
                        },
                    );
                });
            });

            new_matrix
        }

        #[cfg(not(feature = "disk"))]
        {
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
    }
}

impl<T: MatrixElem> Debug for GenericMatrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        #[cfg(feature = "disk")]
        let fmt = f
            .debug_struct("GenericMatrix")
            .field("params", &self.params)
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("file", &self.file)
            .finish();
        #[cfg(not(feature = "disk"))]
        let fmt = f
            .debug_struct("GenericMatrix")
            .field("params", &self.params)
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("inner", &self.inner)
            .finish();
        fmt
    }
}

#[cfg(feature = "disk")]
impl<T: MatrixElem> Clone for GenericMatrix<T> {
    fn clone(&self) -> Self {
        let mut new_matrix = Self::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            self.block_entries(row_offsets, col_offsets)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

#[cfg(feature = "disk")]
impl<T: MatrixElem> PartialEq for GenericMatrix<T> {
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

#[cfg(not(feature = "disk"))]
impl<T: MatrixElem> PartialEq for GenericMatrix<T> {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner &&
            self.params == other.params &&
            self.nrow == other.nrow &&
            self.ncol == other.ncol
    }
}

impl<T: MatrixElem> Eq for GenericMatrix<T> {}

impl<T: MatrixElem> Drop for GenericMatrix<T> {
    fn drop(&mut self) {
        #[cfg(feature = "disk")]
        self.file.set_len(0).expect("failed to truncate file");
    }
}

impl<T: MatrixElem> Add for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<T: MatrixElem> Add<&GenericMatrix<T>> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn add(self, rhs: &Self) -> Self::Output {
        &self + rhs
    }
}

impl<T: MatrixElem> Add<&GenericMatrix<T>> for &GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn add(self, rhs: &GenericMatrix<T>) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );

        let mut new_matrix = GenericMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            add_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MatrixElem> Sub<GenericMatrix<T>> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<T: MatrixElem> Sub<&GenericMatrix<T>> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn sub(self, rhs: &Self) -> Self::Output {
        &self - rhs
    }
}

impl<T: MatrixElem> Sub<&GenericMatrix<T>> for &GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn sub(self, rhs: &GenericMatrix<T>) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );

        let mut new_matrix = GenericMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<T>> {
            let self_block_polys = self.block_entries(row_offsets.clone(), col_offsets.clone());
            let rhs_block_polys = rhs.block_entries(row_offsets, col_offsets);
            sub_block_matrices(self_block_polys, &rhs_block_polys)
        };
        new_matrix.replace_entries(0..self.nrow, 0..self.ncol, f);
        new_matrix
    }
}

impl<T: MatrixElem> Mul<GenericMatrix<T>> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&GenericMatrix<T>> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: &Self) -> Self::Output {
        &self * rhs
    }
}

impl<T: MatrixElem> Mul<GenericMatrix<T>> for &GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: GenericMatrix<T>) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&GenericMatrix<T>> for &GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: &GenericMatrix<T>) -> Self::Output {
        debug_assert!(
            self.ncol == rhs.nrow,
            "Multiplication condition failed: self.ncol ({}) must equal rhs.nrow ({})",
            self.ncol,
            rhs.nrow
        );

        let mut new_matrix = GenericMatrix::new_empty(&self.params, self.nrow, rhs.ncol);
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

impl<T: MatrixElem> Mul<T> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;
    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&T> for GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        &self * rhs
    }
}

impl<T: MatrixElem> Mul<T> for &GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: T) -> Self::Output {
        self * &rhs
    }
}

impl<T: MatrixElem> Mul<&T> for &GenericMatrix<T> {
    type Output = GenericMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        let mut new_matrix = GenericMatrix::new_empty(&self.params, self.nrow, self.ncol);
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

impl<T: MatrixElem> Neg for GenericMatrix<T> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut new_matrix = GenericMatrix::new_empty(&self.params, self.nrow, self.ncol);
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

#[cfg(feature = "disk")]
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

#[cfg(feature = "disk")]
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
