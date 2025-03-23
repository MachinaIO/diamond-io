use super::{DCRTPoly, DCRTPolyParams};
use crate::{
    impl_binop_with_refs, parallel_iter,
    poly::{Poly, PolyMatrix, PolyParams},
};
use itertools::Itertools;
#[cfg(target_os = "linux")]
use memmap2::RemapOptions;
use memmap2::{Advice, Mmap, MmapMut, MmapOptions};
use once_cell::sync::OnceCell;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::{
    fmt::Debug,
    fs::{File, OpenOptions},
    marker::PhantomData,
    ops::{Add, Deref, DerefMut, Index, IndexMut, Mul, Neg, Range, Sub},
    path::{Path, PathBuf},
};
use sysinfo::System;
use tempfile::tempfile;

static BLOCK_SIZE: OnceCell<usize> = OnceCell::new();

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Volatility {
    Persistent,
    Transient,
}

pub struct DCRTPolyMatrix {
    params: DCRTPolyParams,
    file: File,
    nrow: usize,
    ncol: usize,
    volatile: Volatility,
}

impl PolyMatrix for DCRTPolyMatrix {
    type P = DCRTPoly;

    fn from_poly_vec(params: &DCRTPolyParams, vec: Vec<Vec<DCRTPoly>>) -> Self {
        let nrow = vec.len();
        let ncol = vec[0].len();
        let storage = Self::new_empty(params, nrow, ncol);

        let (row_offsets, col_offsets) = block_offsets(0..nrow, 0..ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                let block_row_len = next_block_row_idx - cur_block_row_idx;
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        let mut new_entries = Vec::with_capacity(block_row_len);
                        for i in 0..block_row_len {
                            new_entries.push(Vec::with_capacity(block_col_len));
                            for j in 0..block_col_len {
                                new_entries[i].push(
                                    vec[cur_block_row_idx + i][cur_block_col_idx + j].clone(),
                                );
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            storage.replace_block_entries(
                                *cur_block_col_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_entries,
                            );
                        }
                    },
                );
            },
        );
        storage
    }

    fn entry(&self, i: usize, j: usize) -> Self::P {
        self.block_entries(i..i + 1, j..j + 1)[0][0].clone()
    }

    fn get_row(&self, i: usize) -> Vec<Self::P> {
        self.block_entries(i..i + 1, 0..self.ncol)[0].clone()
    }

    fn get_column(&self, j: usize) -> Vec<Self::P> {
        self.block_entries(0..self.nrow, j..j + 1).iter().map(|row| row[0].clone()).collect()
    }

    fn size(&self) -> (usize, usize) {
        (self.nrow, self.ncol)
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
        let new_matrix = Self::new_empty(&self.params, nrow, ncol);
        let (row_offsets, col_offsets) =
            block_offsets(row_start..row_end, column_start..column_end);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let new_entries = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_entries,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }

    fn zero(params: &<Self::P as Poly>::Params, nrow: usize, ncol: usize) -> Self {
        Self::new_empty(params, nrow, ncol)
    }

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        let new_matrix = Self::new_empty(params, size, size);
        let (row_offsets, col_offsets) = block_offsets(0..size, 0..size);
        let scalar = scalar.unwrap_or_else(|| DCRTPoly::const_one(params));
        let zero_elem = DCRTPoly::const_zero(params);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let mut new_entries =
                            Vec::with_capacity(next_block_row_idx - cur_block_row_idx);
                        for i in 0..(next_block_row_idx - cur_block_row_idx) {
                            new_entries
                                .push(Vec::with_capacity(next_block_col_idx - cur_block_col_idx));
                            for j in 0..(next_block_col_idx - cur_block_col_idx) {
                                if i == j {
                                    new_entries[i].push(scalar.clone());
                                } else {
                                    new_entries[i].push(zero_elem.clone());
                                }
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_entries,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }

    fn transpose(&self) -> Self {
        let new_matrix = Self::new_empty(&self.params, self.ncol, self.nrow);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let cur_entries = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let mut new_entries =
                            Vec::with_capacity(next_block_col_idx - cur_block_col_idx);
                        for j in 0..(next_block_col_idx - cur_block_col_idx) {
                            new_entries
                                .push(Vec::with_capacity(next_block_row_idx - cur_block_row_idx));
                            for i in 0..(next_block_row_idx - cur_block_row_idx) {
                                new_entries[j].push(cur_entries[i][j].clone());
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_col_idx..*next_block_row_idx,
                                *cur_block_row_idx..*next_block_col_idx,
                                new_entries,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
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
        let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);

        let new_matrix = Self::new_empty(&self.params, self.nrow, updated_ncol);
        let (row_offsets, self_col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        let others_col_offsets =
            others.iter().map(|other| block_offsets(0..other.nrow, 0..other.ncol).1).collect_vec();

        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(self_col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );

                let mut col_acc = self.ncol;

                for (other_idx, other) in others.iter().enumerate() {
                    let other_col_offsets = &others_col_offsets[other_idx];
                    parallel_iter!(other_col_offsets.iter().tuple_windows().collect_vec())
                        .for_each(|(cur_block_col_idx, next_block_col_idx)| {
                            let other_block_polys = other.block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                            );
                            // This is secure because the modified entries are not overlapped among
                            // threads
                            unsafe {
                                new_matrix.replace_block_entries(
                                    *cur_block_row_idx..*next_block_row_idx,
                                    col_acc + cur_block_col_idx..col_acc + next_block_col_idx,
                                    other_block_polys,
                                );
                            }
                        });
                    col_acc += other.ncol;
                }
                debug_assert_eq!(col_acc, updated_ncol);
            },
        );
        new_matrix
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[&Self]) -> Self {
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

        let new_matrix = Self::new_empty(&self.params, updated_nrow, self.ncol);
        let (self_row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        let others_row_offsets =
            others.iter().map(|other| block_offsets(0..other.nrow, 0..other.ncol).0).collect_vec();

        // Copy blocks from self
        parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_col_idx, next_block_col_idx)| {
                parallel_iter!(self_row_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_row_idx, next_block_row_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );

                // Copy blocks from others
                let mut row_acc = self.nrow;
                for (other_idx, other) in others.iter().enumerate() {
                    let other_row_offsets = &others_row_offsets[other_idx];
                    parallel_iter!(other_row_offsets.iter().tuple_windows().collect_vec())
                        .for_each(|(cur_block_row_idx, next_block_row_idx)| {
                            let other_block_polys = other.block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                            );
                            // This is secure because the modified entries are not overlapped among
                            // threads
                            unsafe {
                                new_matrix.replace_block_entries(
                                    row_acc + cur_block_row_idx..row_acc + next_block_row_idx,
                                    *cur_block_col_idx..*next_block_col_idx,
                                    other_block_polys,
                                );
                            }
                        });
                    row_acc += other.nrow;
                }
                debug_assert_eq!(row_acc, updated_nrow);
            },
        );

        new_matrix
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[&Self]) -> Self {
        let updated_nrow = others.iter().fold(self.nrow, |acc, other| acc + other.nrow);
        let updated_ncol = others.iter().fold(self.ncol, |acc, other| acc + other.ncol);

        let new_matrix = Self::new_empty(&self.params, updated_nrow, updated_ncol);

        // Copy blocks from self to the top-left corner
        let (self_row_offsets, self_col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(self_row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(self_col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );
            },
        );

        // Copy blocks from others to the diagonal
        let mut row_acc = self.nrow;
        let mut col_acc = self.ncol;

        for other in others {
            let (other_row_offsets, other_col_offsets) =
                block_offsets(0..other.nrow, 0..other.ncol);

            parallel_iter!(other_row_offsets.iter().tuple_windows().collect_vec()).for_each(
                |(cur_block_row_idx, next_block_row_idx)| {
                    parallel_iter!(other_col_offsets.iter().tuple_windows().collect_vec())
                        .for_each(|(cur_block_col_idx, next_block_col_idx)| {
                            let other_block_polys = other.block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                            );
                            // This is secure because the modified entries are not overlapped among
                            // threads
                            unsafe {
                                new_matrix.replace_block_entries(
                                    row_acc + cur_block_row_idx..row_acc + next_block_row_idx,
                                    col_acc + cur_block_col_idx..col_acc + next_block_col_idx,
                                    other_block_polys,
                                );
                            }
                        });
                },
            );

            row_acc += other.nrow;
            col_acc += other.ncol;
        }

        new_matrix
    }

    fn tensor(&self, other: &Self) -> Self {
        let mut new_matrix =
            Self::new_empty(&self.params, self.nrow * other.nrow, self.ncol * other.ncol);
        // for i in 0..self.nrow {
        //     for j in 0..self.ncol {
        //         let scalar = self.entry(i, j);
        //         let sub_matrix = other.clone() * scalar;

        //     }
        // }
        new_matrix
    }

    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self {
        todo!();
    }

    fn decompose(&self) -> Self {
        todo!();
    }

    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self {
        todo!();
    }

    fn mul_tensor_identity(&self, other: &Self, identity_size: usize) -> Self {
        todo!();
    }

    fn mul_tensor_identity_decompose(&self, other: &Self, identity_size: usize) -> Self {
        todo!();
    }

    fn get_column_matrix_decompose(&self, j: usize) -> Self {
        todo!();
    }
}

impl Debug for DCRTPolyMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix")
            .field("params", &self.params)
            .field("nrow", &self.nrow)
            .field("ncol", &self.ncol)
            .field("file", &self.file)
            .finish()
    }
}

impl Clone for DCRTPolyMatrix {
    fn clone(&self) -> Self {
        let new_matrix = Self::new_empty(&self.params, self.nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }
}

impl PartialEq for DCRTPolyMatrix {
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

impl Eq for DCRTPolyMatrix {}

impl DCRTPolyMatrix {
    fn new_empty(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        let entry_size = params.modulus_bits().div_ceil(8);
        let len = entry_size * nrow * ncol;
        let file = tempfile().expect("failed to open file");
        file.set_len(len as u64).expect("failed to set file length");
        if BLOCK_SIZE.get().is_none() {
            let system = System::new();
            let mem_size = system.total_memory() * 2 / 3;
            #[cfg(feature = "parallel")]
            let num_threads = rayon::current_num_threads();
            #[cfg(not(feature = "parallel"))]
            let num_threads = 1;
            let num_polys = mem_size as usize / num_threads / entry_size;
            let block_size = (num_polys as f64).sqrt() as usize;
            BLOCK_SIZE.set(block_size).unwrap();
        }
        Self { params: params.clone(), file, nrow, ncol, volatile: Volatility::Transient }
    }

    fn entry_size(&self) -> usize {
        self.params.modulus_bits().div_ceil(8)
    }

    fn block_entries(
        &self,
        rows: Range<usize>,
        cols: Range<usize>,
    ) -> Vec<Vec<<Self as PolyMatrix>::P>> {
        let entry_size = self.entry_size();
        parallel_iter!(rows)
            .map(|i| {
                let offset = entry_size * (i * self.ncol + cols.start);
                let mmap = map_file(&self.file, offset, entry_size * cols.len());
                let row_vec = mmap.to_vec();
                let row_col_vec = row_vec
                    .chunks(entry_size)
                    .map(|entry| {
                        <<Self as PolyMatrix>::P as Poly>::from_bytes(&self.params, &entry)
                    })
                    .collect_vec();
                drop(mmap);
                row_col_vec
            })
            .collect::<Vec<Vec<_>>>()
    }

    unsafe fn replace_block_entries(
        &self,
        rows: Range<usize>,
        cols: Range<usize>,
        new_entries: Vec<Vec<<Self as PolyMatrix>::P>>,
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
                .map(|poly| poly.to_bytes())
                .flatten()
                .collect::<Vec<_>>();
            mmap.copy_from_slice(&bytes);
            drop(mmap);
        });
    }

    // // (m * n1), (m * n2) -> (m * (n1 + n2))
    // pub fn concat_columns(&self, other: &Self) -> Self {
    //     debug_assert_eq!(self.nrow, other.nrow);
    //     debug_assert_eq!(self.entry_size, other.entry_size);
    //     let updated_ncol = self.ncol + other.ncol;

    //     let new_matrix = Self::new_empty(self.nrow, updated_ncol, self.entry_size);
    //     parallel_iter!(0..self.nrow).for_each(|i| {
    //         let self_mmap = map_file(
    //             &self.file,
    //             self.entry_size * (i * self.ncol),
    //             self.entry_size * self.ncol,
    //         );
    //         let other_mmap = map_file(
    //             &other.file,
    //             other.entry_size * (i * other.ncol),
    //             other.entry_size * other.ncol,
    //         );
    //         let mut new_mmap = map_file_mut(
    //             &new_matrix.file,
    //             self.entry_size * (i * updated_ncol),
    //             self.entry_size * updated_ncol,
    //         );
    //         new_mmap[0..self.entry_size * self.ncol].copy_from_slice(&self_mmap);
    //         drop(self_mmap);
    //         new_mmap[self.entry_size * self.ncol..].copy_from_slice(&other_mmap);
    //         drop(other_mmap);
    //         drop(new_mmap);
    //     });

    //     new_matrix
    // }

    // // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    // pub fn concat_rows(&self, other: &Self) -> Self {
    //     debug_assert_eq!(self.ncol, other.ncol);
    //     debug_assert_eq!(self.entry_size, other.entry_size);
    //     let updated_nrow = self.nrow + other.nrow;

    //     let new_matrix = Self::new_empty(updated_nrow, self.ncol, self.entry_size);
    //     parallel_iter!(0..self.nrow).for_each(|i| {
    //         let offset = self.entry_size * i * self.ncol;
    //         let self_mmap = map_file(&self.file, offset, self.entry_size * self.ncol);
    //         let mut new_mmap = map_file_mut(&new_matrix.file, offset, self.entry_size *
    // self.ncol);         new_mmap.copy_from_slice(&self_mmap);
    //         drop(self_mmap);
    //         drop(new_mmap);
    //     });
    //     parallel_iter!(0..other.nrow).for_each(|i| {
    //         let other_mmap = map_file(
    //             &other.file,
    //             other.entry_size * i * other.ncol,
    //             other.entry_size * other.ncol,
    //         );
    //         let mut new_mmap = map_file_mut(
    //             &new_matrix.file,
    //             self.entry_size * (self.nrow + i) * self.ncol,
    //             self.entry_size * self.ncol,
    //         );
    //         new_mmap.copy_from_slice(&other_mmap);
    //         drop(other_mmap);
    //         drop(new_mmap);
    //     });
    //     new_matrix
    // }

    // // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    // pub fn concat_diag(&self, other: &Self) -> Self {
    //     debug_assert_eq!(self.entry_size, other.entry_size);
    //     let updated_nrow = self.nrow + other.nrow;
    //     let updated_ncol = self.ncol + other.ncol;

    //     let new_matrix = Self::new_empty(updated_nrow, updated_ncol, self.entry_size);

    //     // Copy the first matrix (self) to the top-left corner
    //     parallel_iter!(0..self.nrow).for_each(|i| {
    //         let self_offset = self.entry_size * i * self.ncol;
    //         let new_offset = self.entry_size * i * updated_ncol;
    //         let self_mmap = map_file(&self.file, self_offset, self.entry_size * self.ncol);
    //         let mut new_mmap =
    //             map_file_mut(&new_matrix.file, new_offset, self.entry_size * self.ncol);
    //         new_mmap.copy_from_slice(&self_mmap);
    //         drop(self_mmap);
    //         drop(new_mmap);
    //     });

    //     // Copy the second matrix (other) to the bottom-right corner
    //     parallel_iter!(0..other.nrow).for_each(|i| {
    //         let other_offset = other.entry_size * i * other.ncol;
    //         let new_offset = self.entry_size * ((self.nrow + i) * updated_ncol + self.ncol);
    //         let other_mmap = map_file(&other.file, other_offset, other.entry_size * other.ncol);
    //         let mut new_mmap =
    //             map_file_mut(&new_matrix.file, new_offset, other.entry_size * other.ncol);
    //         new_mmap.copy_from_slice(&other_mmap);
    //         drop(other_mmap);
    //         drop(new_mmap);
    //     });

    //     new_matrix
    // }
}

impl Drop for DCRTPolyMatrix {
    fn drop(&mut self) {
        if self.volatile == Volatility::Transient {
            self.file.set_len(0).expect("failed to truncate file");
        }
    }
}

impl Add for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl Add<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn add(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );
        let new_matrix = DCRTPolyMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let mut self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let rhs_block_polys = rhs.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let block_row_len = next_block_row_idx - cur_block_row_idx;
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        for i in 0..block_row_len {
                            for j in 0..block_col_len {
                                self_block_polys[i][j] += &rhs_block_polys[i][j];
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }
}

impl Sub for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl Sub<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn sub(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        debug_assert!(
            self.nrow == rhs.nrow && self.ncol == rhs.ncol,
            "Addition requires matrices of same dimensions: self({}, {}) != rhs({}, {})",
            self.nrow,
            self.ncol,
            rhs.nrow,
            rhs.ncol
        );
        let new_matrix = DCRTPolyMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let mut self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let rhs_block_polys = rhs.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let block_row_len = next_block_row_idx - cur_block_row_idx;
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        for i in 0..block_row_len {
                            for j in 0..block_col_len {
                                self_block_polys[i][j] -= &rhs_block_polys[i][j];
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }
}

impl Mul for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&DCRTPolyMatrix> for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn mul(self, rhs: &DCRTPolyMatrix) -> Self::Output {
        debug_assert!(
            self.ncol == rhs.nrow,
            "Multiplication condition failed: self.ncol ({}) must equal rhs.nrow ({})",
            self.ncol,
            rhs.nrow
        );
        let new_matrix = DCRTPolyMatrix::new_empty(&self.params, self.nrow, rhs.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..rhs.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        for ip_block_idx in 0..self.ncol {
                            let mut ip_sum = new_matrix.block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                            );
                            let self_block_polys = self.block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                ip_block_idx..ip_block_idx + 1,
                            );
                            let other_block_polys = rhs.block_entries(
                                ip_block_idx..ip_block_idx + 1,
                                *cur_block_col_idx..*next_block_col_idx,
                            );
                            for i in 0..(next_block_row_idx - cur_block_row_idx) {
                                for j in 0..(next_block_col_idx - cur_block_col_idx) {
                                    let mut sum = DCRTPoly::const_zero(&self.params);
                                    for k in 0..(self.ncol) {
                                        sum += &self_block_polys[i][0] * &other_block_polys[0][j];
                                    }
                                    ip_sum[i][j] += sum;
                                }
                            }
                            // This is secure because the modified entries are not overlapped among
                            // threads
                            unsafe {
                                new_matrix.replace_block_entries(
                                    *cur_block_row_idx..*next_block_row_idx,
                                    *cur_block_col_idx..*next_block_col_idx,
                                    ip_sum,
                                );
                            }
                        }
                    },
                );
            },
        );
        new_matrix
    }
}

impl Mul<DCRTPoly> for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;
    fn mul(self, rhs: DCRTPoly) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&DCRTPoly> for DCRTPolyMatrix {
    type Output = DCRTPolyMatrix;

    fn mul(self, rhs: &DCRTPoly) -> Self::Output {
        let new_matrix = DCRTPolyMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let mut self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let block_row_len = next_block_row_idx - cur_block_row_idx;
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        for i in 0..block_row_len {
                            for j in 0..block_col_len {
                                self_block_polys[i][j] *= rhs;
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }
}

impl Neg for DCRTPolyMatrix {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let new_matrix = DCRTPolyMatrix::new_empty(&self.params, self.nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let mut self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let block_row_len = next_block_row_idx - cur_block_row_idx;
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        for i in 0..block_row_len {
                            for j in 0..block_col_len {
                                self_block_polys[i][j] = -self_block_polys[i][j].clone();
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                self_block_polys,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }
}

// fn dim1_to_dim2(one_dim: usize, ncol: usize) -> (usize, usize) {
//     (one_dim / ncol, one_dim % ncol)
// }

fn map_file(file: &File, offset: usize, len: usize) -> Mmap {
    unsafe {
        MmapOptions::new().offset(offset as u64).len(len).map(file).expect("failed to map file")
    }
}

unsafe fn map_file_mut(file: &File, offset: usize, len: usize) -> MmapMut {
    unsafe {
        MmapOptions::new().offset(offset as u64).len(len).map_mut(file).expect("failed to map file")
    }
}

fn block_offsets(rows: Range<usize>, cols: Range<usize>) -> (Vec<usize>, Vec<usize>) {
    let block_size = *BLOCK_SIZE.get().unwrap();
    let nrow = rows.end - rows.start;
    let ncol = cols.end - cols.start;
    let num_blocks_row = rows.len().div_ceil(block_size);
    let num_blocks_col = cols.len().div_ceil(block_size);
    let row_offsets = (0..num_blocks_row)
        .map(|i| if i == num_blocks_row - 1 { nrow - i * block_size } else { block_size })
        .scan(0, |state, len| {
            let offset = rows.start + *state;
            *state += len;
            Some(offset)
        })
        .collect();
    let col_offsets = (0..num_blocks_col)
        .map(|i| if i == num_blocks_col - 1 { ncol - i * block_size } else { block_size })
        .scan(0, |state, len| {
            let offset = cols.start + *state;
            *state += len;
            Some(offset)
        })
        .collect();
    (row_offsets, col_offsets)
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     // Helper function to create a test matrix with specified dimensions
//     fn create_test_matrix(nrow: usize, ncol: usize) -> MatrixStorage {
//         let mut entries = Vec::with_capacity(nrow);
//         for i in 0..nrow {
//             let mut row = Vec::with_capacity(ncol);
//             for j in 0..ncol {
//                 // Create a unique value for each cell based on its position
//                 let value = ((i * ncol + j) as u8) % 255;
//                 row.push(vec![value, value + 1, value + 2, value + 3]);
//             }
//             entries.push(row);
//         }
//         MatrixStorage::new(entries)
//     }

//     #[test]
//     fn test_matrix_storage_new_and_entries() {
//         // Create a 3x4 matrix
//         let nrow = 3;
//         let ncol = 4;
//         let entry_size = 4; // Each entry is 4 bytes

//         // Create test data
//         let mut entries = Vec::with_capacity(nrow);
//         for i in 0..nrow {
//             let mut row = Vec::with_capacity(ncol);
//             for j in 0..ncol {
//                 let value = ((i * ncol + j) as u8) % 255;
//                 row.push(vec![value, value + 1, value + 2, value + 3]);
//             }
//             entries.push(row);
//         }

//         // Create matrix
//         let matrix = MatrixStorage::new(entries.clone());

//         // Verify dimensions
//         assert_eq!(matrix.nrow, nrow);
//         assert_eq!(matrix.ncol, ncol);
//         assert_eq!(matrix.entry_size, entry_size);

//         // Retrieve all entries and verify they match the original data
//         let retrieved_entries = matrix.entries(0..nrow, 0..ncol);
//         assert_eq!(retrieved_entries, entries);
//     }

//     #[test]
//     fn test_matrix_storage_replace_entries() {
//         // Create a 3x4 matrix
//         let matrix = create_test_matrix(3, 4);

//         // Create new entries for a 2x2 submatrix
//         let new_entries = vec![
//             vec![vec![100, 101, 102, 103], vec![110, 111, 112, 113]],
//             vec![vec![120, 121, 122, 123], vec![130, 131, 132, 133]],
//         ];

//         // Replace a 2x2 submatrix starting at position (1, 1)
//         let mut modified_matrix = matrix.clone();
//         modified_matrix.replace_entries(1..3, 1..3, new_entries.clone());

//         // Retrieve the modified submatrix and verify it matches the new entries
//         let retrieved_entries = modified_matrix.entries(1..3, 1..3);
//         assert_eq!(retrieved_entries, new_entries);

//         // Verify that the rest of the matrix is unchanged
//         let top_left = modified_matrix.entries(0..1, 0..1);
//         assert_eq!(top_left, matrix.entries(0..1, 0..1));

//         let top_right = modified_matrix.entries(0..1, 3..4);
//         assert_eq!(top_right, matrix.entries(0..1, 3..4));

//         let bottom_left = modified_matrix.entries(2..3, 0..1);
//         assert_eq!(bottom_left, matrix.entries(2..3, 0..1));
//     }

//     #[test]
//     fn test_matrix_storage_concat_columns() {
//         // Create two matrices with the same number of rows
//         let matrix1 = create_test_matrix(3, 2);
//         let matrix2 = create_test_matrix(3, 3);

//         // Concatenate the matrices horizontally
//         let result = matrix1.concat_columns(&matrix2);

//         // Verify dimensions
//         assert_eq!(result.nrow, 3);
//         assert_eq!(result.ncol, 5);
//         assert_eq!(result.entry_size, 4);

//         // Verify the left part matches matrix1
//         let left_part = result.entries(0..3, 0..2);
//         assert_eq!(left_part, matrix1.entries(0..3, 0..2));

//         // Verify the right part matches matrix2
//         let right_part = result.entries(0..3, 2..5);
//         assert_eq!(right_part, matrix2.entries(0..3, 0..3));
//     }

//     #[test]
//     fn test_matrix_storage_concat_rows() {
//         // Create two matrices with the same number of columns
//         let matrix1 = create_test_matrix(2, 3);
//         let matrix2 = create_test_matrix(3, 3);

//         // Concatenate the matrices vertically
//         let result = matrix1.concat_rows(&matrix2);

//         // Verify dimensions
//         assert_eq!(result.nrow, 5);
//         assert_eq!(result.ncol, 3);
//         assert_eq!(result.entry_size, 4);

//         // Verify the top part matches matrix1
//         let top_part = result.entries(0..2, 0..3);
//         assert_eq!(top_part, matrix1.entries(0..2, 0..3));

//         // Verify the bottom part matches matrix2
//         let bottom_part = result.entries(2..5, 0..3);
//         assert_eq!(bottom_part, matrix2.entries(0..3, 0..3));
//     }

//     #[test]
//     fn test_matrix_storage_concat_diag() {
//         // Create two matrices
//         let matrix1 = create_test_matrix(2, 3);
//         let matrix2 = create_test_matrix(3, 2);

//         // Concatenate the matrices diagonally
//         let result = matrix1.concat_diag(&matrix2);

//         // Verify dimensions
//         assert_eq!(result.nrow, 5);
//         assert_eq!(result.ncol, 5);
//         assert_eq!(result.entry_size, 4);

//         // Verify the top-left part matches matrix1
//         let top_left = result.entries(0..2, 0..3);
//         assert_eq!(top_left, matrix1.entries(0..2, 0..3));

//         // Verify the bottom-right part matches matrix2
//         let bottom_right = result.entries(2..5, 3..5);
//         assert_eq!(bottom_right, matrix2.entries(0..3, 0..2));

//         // Verify the top-right part is zeros (or empty)
//         let top_right = result.entries(0..2, 3..5);
//         for row in top_right {
//             for entry in row {
//                 assert_eq!(entry, vec![0, 0, 0, 0]);
//             }
//         }

//         // Verify the bottom-left part is zeros (or empty)
//         let bottom_left = result.entries(2..5, 0..3);
//         for row in bottom_left {
//             for entry in row {
//                 assert_eq!(entry, vec![0, 0, 0, 0]);
//             }
//         }
//     }

//     // #[test]
//     // fn test_matrix_storage_volatility() {
//     //     // Create a matrix with default transient volatility
//     //     let mut matrix = create_test_matrix(2, 2);
//     //     assert_eq!(matrix.volatile, Volatility::Transient);

//     //     // Change to persistent
//     //     matrix.to_persistent();
//     //     assert_eq!(matrix.volatile, Volatility::Persistent);

//     //     // Change back to transient
//     //     matrix.to_transient();
//     //     assert_eq!(matrix.volatile, Volatility::Transient);
//     // }
// }
