use super::{DCRTPoly, DCRTPolyParams, FinRingElem};
use crate::{
    parallel_iter,
    poly::{Poly, PolyMatrix, PolyParams},
};
use itertools::Itertools;
#[cfg(target_os = "linux")]
use memmap2::RemapOptions;
use memmap2::{Mmap, MmapMut, MmapOptions};
use num_bigint::BigInt;
use once_cell::sync::OnceCell;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::{
    fmt::Debug,
    fs::File,
    ops::{Add, Mul, Neg, Range, Sub},
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
                                *cur_block_row_idx - row_start..*next_block_row_idx - row_start,
                                *cur_block_col_idx - column_start
                                    ..*next_block_col_idx - column_start,
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
        debug_assert_eq!(row_offsets, col_offsets);
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
        let new_matrix =
            Self::new_empty(&self.params, self.nrow * other.nrow, self.ncol * other.ncol);
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                let scalar = self.entry(i, j);
                let sub_matrix = other.clone() * scalar;
                let (row_offsets, col_offsets) =
                    block_offsets(0..sub_matrix.nrow, 0..sub_matrix.ncol);
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
                                        i * sub_matrix.nrow + *cur_block_row_idx
                                            ..i * sub_matrix.nrow + *next_block_row_idx,
                                        j * sub_matrix.ncol + *cur_block_col_idx
                                            ..j * sub_matrix.ncol + *next_block_col_idx,
                                        sub_block_polys,
                                    );
                                }
                            },
                        );
                    },
                );
            }
        }
        new_matrix
    }

    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self {
        let bit_length = params.modulus_bits();
        let mut poly_vec = Vec::with_capacity(bit_length);
        for i in 0u32..(bit_length as u32) {
            let value = BigInt::from(2).pow(i);
            poly_vec.push(DCRTPoly::from_const(params, &FinRingElem::new(value, params.modulus())));
        }
        let gadget_vector = Self::from_poly_vec(params, vec![poly_vec]);
        let identity = DCRTPolyMatrix::identity(params, size, None);
        identity.tensor(&gadget_vector)
    }

    fn decompose(&self) -> Self {
        let bit_length = self.params.modulus_bits();
        let new_nrow = self.nrow * bit_length;
        let new_matrix = Self::new_empty(&self.params, new_nrow, self.ncol);
        let (row_offsets, col_offsets) = block_offsets(0..self.nrow, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let block_row_len = next_block_row_idx - cur_block_row_idx;
                        let block_col_len = next_block_col_idx - cur_block_col_idx;
                        let mut new_entries =
                            vec![
                                vec![DCRTPoly::const_zero(&self.params); block_col_len];
                                block_row_len * bit_length
                            ];
                        for j in 0..block_col_len {
                            for i in 0..block_row_len {
                                let poly = &self_block_polys[i][j];
                                let decomposed_polys = poly.decompose(&self.params);
                                for (k, poly) in decomposed_polys.into_iter().enumerate() {
                                    new_entries[i * bit_length + k][j] = poly;
                                }
                            }
                        }
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                cur_block_row_idx * bit_length..next_block_row_idx * bit_length,
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

    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self {
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
                        let new_block_polys = self_block_polys
                            .iter()
                            .map(|row| {
                                row.iter()
                                    .map(|poly| {
                                        poly.modulus_switch(&self.params, new_modulus.clone())
                                    })
                                    .collect_vec()
                            })
                            .collect_vec();
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_block_polys,
                            );
                        }
                    },
                );
            },
        );
        new_matrix
    }

    fn mul_tensor_identity(&self, other: &Self, identity_size: usize) -> Self {
        debug_assert_eq!(self.ncol, other.nrow * identity_size);

        let slice_width = other.nrow;
        let mut slice_results = Vec::with_capacity(identity_size);

        for i in 0..identity_size {
            let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
            slice_results.push(slice * other);
        }

        slice_results[0].clone().concat_columns(&slice_results[1..].iter().collect::<Vec<_>>())
    }

    fn mul_tensor_identity_decompose(&self, other: &Self, identity_size: usize) -> Self {
        debug_assert_eq!(self.ncol, other.nrow * identity_size * self.params.modulus_bits());

        let slice_width = other.nrow * self.params.modulus_bits();
        let mut output = Vec::with_capacity(other.ncol * identity_size);

        for j in 0..other.ncol {
            let jth_col_m_decompose = other.get_column_matrix_decompose(j);
            for i in 0..identity_size {
                let slice = self.slice(0, self.nrow, i * slice_width, (i + 1) * slice_width);
                let result = slice * &jth_col_m_decompose;
                output.push(result);
            }
        }

        output[0].clone().concat_columns(&output[1..].iter().collect::<Vec<_>>())
    }

    fn get_column_matrix_decompose(&self, j: usize) -> Self {
        let column = self.get_column(j);
        let column_matrix =
            Self::from_poly_vec(&self.params, column.into_iter().map(|poly| vec![poly]).collect());
        column_matrix.decompose()
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
                        let self_block_polys = self.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let rhs_block_polys = rhs.block_entries(
                            *cur_block_row_idx..*next_block_row_idx,
                            *cur_block_col_idx..*next_block_col_idx,
                        );
                        let new_block_polys = add_block_matrices(self_block_polys, rhs_block_polys);
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_block_polys,
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
                        let new_block_polys = sub_block_matrices(self_block_polys, rhs_block_polys);
                        // This is secure because the modified entries are not overlapped among
                        // threads
                        unsafe {
                            new_matrix.replace_block_entries(
                                *cur_block_row_idx..*next_block_row_idx,
                                *cur_block_col_idx..*next_block_col_idx,
                                new_block_polys,
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
        let (_, ip_offsets) = block_offsets(0..0, 0..self.ncol);
        parallel_iter!(row_offsets.iter().tuple_windows().collect_vec()).for_each(
            |(cur_block_row_idx, next_block_row_idx)| {
                parallel_iter!(col_offsets.iter().tuple_windows().collect_vec()).for_each(
                    |(cur_block_col_idx, next_block_col_idx)| {
                        // parallel_iter!()
                        parallel_iter!(ip_offsets.iter().tuple_windows().collect_vec()).for_each(
                            |(cur_block_ip_idx, next_block_ip_idx)| {
                                let ip_sum = new_matrix.block_entries(
                                    *cur_block_row_idx..*next_block_row_idx,
                                    *cur_block_col_idx..*next_block_col_idx,
                                );
                                let self_block_polys = self.block_entries(
                                    *cur_block_row_idx..*next_block_row_idx,
                                    *cur_block_ip_idx..*next_block_ip_idx,
                                );
                                let other_block_polys = rhs.block_entries(
                                    *cur_block_ip_idx..*next_block_ip_idx,
                                    *cur_block_col_idx..*next_block_col_idx,
                                );
                                let muled = mul_block_matrices(
                                    &self.params,
                                    self_block_polys,
                                    other_block_polys,
                                );
                                let added = add_block_matrices(ip_sum, muled);
                                // This is secure because the modified entries are not overlapped among
                                // threads
                                unsafe {
                                    new_matrix.replace_block_entries(
                                        *cur_block_row_idx..*next_block_row_idx,
                                        *cur_block_col_idx..*next_block_col_idx,
                                        added,
                                    );
                                }
                            },
                        );
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

impl DCRTPolyMatrix {
    fn new_empty(params: &DCRTPolyParams, nrow: usize, ncol: usize) -> Self {
        let dim = params.ring_dimension() as usize;
        let log_q_bytes = params.modulus_bits().div_ceil(8);
        let entry_size = dim * log_q_bytes;
        let len = entry_size * nrow * ncol;
        let file = tempfile().expect("failed to open file");
        file.set_len(len as u64).expect("failed to set file length");
        if BLOCK_SIZE.get().is_none() {
            let system = System::new_all();
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
        let log_q_bytes = self.params.modulus_bits().div_ceil(8);
        let dim = self.params.ring_dimension() as usize;
        dim * log_q_bytes
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
}

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

fn add_block_matrices(lhs: Vec<Vec<DCRTPoly>>, rhs: Vec<Vec<DCRTPoly>>) -> Vec<Vec<DCRTPoly>> {
    let nrow = lhs.len();
    let ncol = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol)
                .map(|j| lhs[i][j].clone() + &rhs[i][j])
                .collect::<Vec<DCRTPoly>>()
        })
        .collect::<Vec<Vec<DCRTPoly>>>()
}

fn sub_block_matrices(lhs: Vec<Vec<DCRTPoly>>, rhs: Vec<Vec<DCRTPoly>>) -> Vec<Vec<DCRTPoly>> {
    let nrow = lhs.len();
    let ncol = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol)
                .map(|j| lhs[i][j].clone() - &rhs[i][j])
                .collect::<Vec<DCRTPoly>>()
        })
        .collect::<Vec<Vec<DCRTPoly>>>()
}

fn mul_block_matrices(
    params: &DCRTPolyParams,
    lhs: Vec<Vec<DCRTPoly>>,
    rhs: Vec<Vec<DCRTPoly>>,
) -> Vec<Vec<DCRTPoly>> {
    let nrow = lhs.len();
    let ncol = rhs[0].len();
    let n_inner = lhs[0].len();
    parallel_iter!(0..nrow)
        .map(|i| {
            parallel_iter!(0..ncol)
                .map(|j| {
                    let mut sum = DCRTPoly::const_zero(params);
                    for k in 0..n_inner {
                        sum += &lhs[i][k] * &rhs[k][j];
                    }
                    sum
                })
                .collect::<Vec<DCRTPoly>>()
        })
        .collect::<Vec<Vec<DCRTPoly>>>()
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use std::sync::Arc;

    use num_bigint::BigUint;

    use super::*;
    use crate::poly::{
        dcrt::{DCRTPolyParams, DCRTPolyUniformSampler},
        sampler::PolyUniformSampler,
    };

    #[test]
    fn test_matrix_gadget_matrix() {
        let params = DCRTPolyParams::default();
        let size = 3;
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, size);
        assert_eq!(gadget_matrix.size().0, size);
        assert_eq!(gadget_matrix.size().1, size * params.modulus_bits());
    }

    #[test]
    fn test_matrix_decompose() {
        let params = DCRTPolyParams::default();
        let bit_length = params.modulus_bits();

        // Create a simple 2x8 matrix with some non-zero values
        let mut matrix_vec = Vec::with_capacity(2);
        let value = FinRingElem::new(5u32, params.modulus());

        // Create first row
        let mut row1 = Vec::with_capacity(8);
        row1.push(DCRTPoly::from_const(&params, &value));
        for _ in 1..8 {
            row1.push(DCRTPoly::const_zero(&params));
        }

        // Create second row
        let mut row2 = Vec::with_capacity(8);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::from_const(&params, &value));
        for _ in 2..8 {
            row2.push(DCRTPoly::const_zero(&params));
        }

        matrix_vec.push(row1);
        matrix_vec.push(row2);

        let matrix = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
        assert_eq!(matrix.size().0, 2);
        assert_eq!(matrix.size().1, 8);

        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, 2);
        assert_eq!(gadget_matrix.size().0, 2);
        assert_eq!(gadget_matrix.size().1, 2 * bit_length);

        let decomposed = matrix.decompose();
        assert_eq!(decomposed.size().0, 2 * bit_length);
        assert_eq!(decomposed.size().1, 8);

        let expected_matrix = gadget_matrix * decomposed;
        assert_eq!(expected_matrix.size().0, 2);
        assert_eq!(expected_matrix.size().1, 8);
        assert_eq!(matrix, expected_matrix);
    }

    #[test]
    fn test_matrix_basic_operations() {
        let params = DCRTPolyParams::default();

        // Test zero and identity matrices
        let zero = DCRTPolyMatrix::zero(&params, 2, 2);
        let identity = DCRTPolyMatrix::identity(&params, 2, None);

        // Test matrix creation and equality
        let value = FinRingElem::new(5u32, params.modulus());

        // Create a 2x2 matrix with values at (0,0) and (1,1)
        let mut matrix_vec = Vec::with_capacity(2);
        let mut row1 = Vec::with_capacity(2);
        row1.push(DCRTPoly::from_const(&params, &value));
        row1.push(DCRTPoly::const_zero(&params));

        let mut row2 = Vec::with_capacity(2);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::from_const(&params, &value));

        matrix_vec.push(row1);
        matrix_vec.push(row2);

        let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
        assert_eq!(matrix1.entry(0, 0).coeffs()[0], value);
        let matrix2 = matrix1.clone();
        assert_eq!(matrix1, matrix2);

        // Test addition
        let sum = matrix1.clone() + &matrix2;
        let value_10 = FinRingElem::new(10u32, params.modulus());
        assert_eq!(sum.entry(0, 0).coeffs()[0], value_10);

        // Test subtraction
        let diff = matrix1.clone() - &matrix2;
        assert_eq!(diff, zero);

        // Test multiplication
        let prod = matrix1 * &identity;
        assert_eq!(prod.size(), (2, 2));
        // Check that the product has the same values as the original matrix
        assert_eq!(prod.entry(0, 0).coeffs()[0], value);
        assert_eq!(prod.entry(1, 1).coeffs()[0], value);
    }

    #[test]
    fn test_matrix_concatenation() {
        let params = DCRTPolyParams::default();
        let value = FinRingElem::new(5u32, params.modulus());

        // Create first matrix with value at (0,0)
        let mut matrix1_vec = Vec::with_capacity(2);
        let mut row1 = Vec::with_capacity(2);
        row1.push(DCRTPoly::from_const(&params, &value));
        row1.push(DCRTPoly::const_zero(&params));

        let mut row2 = Vec::with_capacity(2);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::const_zero(&params));

        matrix1_vec.push(row1);
        matrix1_vec.push(row2);

        let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix1_vec);

        // Create second matrix with value at (1,1)
        let mut matrix2_vec = Vec::with_capacity(2);
        let mut row1 = Vec::with_capacity(2);
        row1.push(DCRTPoly::const_zero(&params));
        row1.push(DCRTPoly::const_zero(&params));

        let mut row2 = Vec::with_capacity(2);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::from_const(&params, &value));

        matrix2_vec.push(row1);
        matrix2_vec.push(row2);

        let matrix2 = DCRTPolyMatrix::from_poly_vec(&params, matrix2_vec);

        // Test column concatenation
        let col_concat = matrix1.concat_columns(&[&matrix2]);
        assert_eq!(col_concat.size().0, 2);
        assert_eq!(col_concat.size().1, 4);
        assert_eq!(col_concat.entry(0, 0).coeffs()[0], value);
        assert_eq!(col_concat.entry(1, 3).coeffs()[0], value);

        // Test row concatenation
        let row_concat = matrix1.concat_rows(&[&matrix2]);
        assert_eq!(row_concat.size().0, 4);
        assert_eq!(row_concat.size().1, 2);
        assert_eq!(row_concat.entry(0, 0).coeffs()[0], value);
        assert_eq!(row_concat.entry(3, 1).coeffs()[0], value);

        // Test diagonal concatenation
        let diag_concat = matrix1.concat_diag(&[&matrix2]);
        assert_eq!(diag_concat.size().0, 4);
        assert_eq!(diag_concat.size().1, 4);
        assert_eq!(diag_concat.entry(0, 0).coeffs()[0], value);
        assert_eq!(diag_concat.entry(3, 3).coeffs()[0], value);
    }

    #[test]
    fn test_matrix_tensor_product() {
        let params = DCRTPolyParams::default();
        let value = FinRingElem::new(5u32, params.modulus());

        // Create first matrix with value at (0,0)
        let mut matrix1_vec = Vec::with_capacity(2);
        let mut row1 = Vec::with_capacity(2);
        row1.push(DCRTPoly::from_const(&params, &value));
        row1.push(DCRTPoly::const_zero(&params));

        let mut row2 = Vec::with_capacity(2);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::const_zero(&params));

        matrix1_vec.push(row1);
        matrix1_vec.push(row2);

        let matrix1 = DCRTPolyMatrix::from_poly_vec(&params, matrix1_vec);

        // Create second matrix with value at (0,0)
        let mut matrix2_vec = Vec::with_capacity(2);
        let mut row1 = Vec::with_capacity(2);
        row1.push(DCRTPoly::from_const(&params, &value));
        row1.push(DCRTPoly::const_zero(&params));

        let mut row2 = Vec::with_capacity(2);
        row2.push(DCRTPoly::const_zero(&params));
        row2.push(DCRTPoly::const_zero(&params));

        matrix2_vec.push(row1);
        matrix2_vec.push(row2);

        let matrix2 = DCRTPolyMatrix::from_poly_vec(&params, matrix2_vec);

        let tensor = matrix1.tensor(&matrix2);
        assert_eq!(tensor.size().0, 4);
        assert_eq!(tensor.size().1, 4);

        // Check that the (0,0) element is the product of the (0,0) elements
        let value_25 = FinRingElem::new(25u32, params.modulus());
        assert_eq!(tensor.entry(0, 0).coeffs()[0], value_25);
    }

    #[test]
    fn test_matrix_modulus_switch() {
        let params = DCRTPolyParams::default();

        let value00 = FinRingElem::new(1023782870921908217643761278891282178u128, params.modulus());
        let value01 = FinRingElem::new(8179012198875468938912873783289218738u128, params.modulus());
        let value10 = FinRingElem::new(2034903202902173762872163465127672178u128, params.modulus());
        let value11 = FinRingElem::new(1990091289902891278121564387120912660u128, params.modulus());

        let matrix_vec = vec![
            vec![DCRTPoly::from_const(&params, &value00), DCRTPoly::from_const(&params, &value01)],
            vec![DCRTPoly::from_const(&params, &value10), DCRTPoly::from_const(&params, &value11)],
        ];

        let matrix = DCRTPolyMatrix::from_poly_vec(&params, matrix_vec);
        let new_modulus = Arc::new(BigUint::from(2u32));
        let switched = matrix.modulus_switch(&new_modulus);

        // Although the value becomes less than the new modulus, the set modulus is still the same
        assert_eq!(switched.params.modulus(), params.modulus());

        let new_value00 = value00.modulus_switch(new_modulus.clone());
        let new_value01 = value01.modulus_switch(new_modulus.clone());
        let new_value10 = value10.modulus_switch(new_modulus.clone());
        let new_value11 = value11.modulus_switch(new_modulus.clone());

        let expected_vec = vec![
            vec![
                DCRTPoly::from_const(&params, &new_value00),
                DCRTPoly::from_const(&params, &new_value01),
            ],
            vec![
                DCRTPoly::from_const(&params, &new_value10),
                DCRTPoly::from_const(&params, &new_value11),
            ],
        ];

        let expected = DCRTPolyMatrix::from_poly_vec(&params, expected_vec);
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
        assert_eq!(result.size().0, 2);
        assert_eq!(result.size().1, 12);

        let identity = DCRTPolyMatrix::identity(&params, 4, None);

        // Check result
        let expected_result = s * (identity.tensor(&other));

        assert_eq!(expected_result.size().0, 2);
        assert_eq!(expected_result.size().1, 12);
        assert_eq!(result, expected_result)
    }

    #[test]
    fn test_matrix_mul_tensor_identity_decompose_naive() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix S (2x2516)
        let s =
            sampler.sample_uniform(&params, 2, 2516, crate::poly::sampler::DistType::FinRingDist);

        // Create 'other' matrix (2x68)
        let other =
            sampler.sample_uniform(&params, 2, 68, crate::poly::sampler::DistType::FinRingDist);

        // Decompose 'other' matrix
        let other_decompose = other.decompose();

        // Perform S * (I_37 ⊗ G^-1(other))
        let result: DCRTPolyMatrix = s.mul_tensor_identity(&other_decompose, 37);

        // Check dimensions
        assert_eq!(result.size().0, 2);
        assert_eq!(result.size().1, 2516);

        let identity = DCRTPolyMatrix::identity(&params, 37, None);

        // Check result
        let expected_result = s * (identity.tensor(&other_decompose));

        assert_eq!(expected_result.size().0, 2);
        assert_eq!(expected_result.size().1, 2516);
        assert_eq!(result, expected_result)
    }

    #[test]
    fn test_matrix_mul_tensor_identity_decompose_optimal() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Create matrix S (2x2516)
        let s =
            sampler.sample_uniform(&params, 2, 2516, crate::poly::sampler::DistType::FinRingDist);

        // Create 'other' matrix (2x68)
        let other =
            sampler.sample_uniform(&params, 2, 68, crate::poly::sampler::DistType::FinRingDist);

        // Perform S * (I_37 ⊗ G^-1(other))
        let result: DCRTPolyMatrix = s.mul_tensor_identity_decompose(&other, 37);

        // Check dimensions
        assert_eq!(result.size().0, 2);
        assert_eq!(result.size().1, 2516);

        let identity = DCRTPolyMatrix::identity(&params, 37, None);

        // Check result
        let expected_result = s.clone() * (identity.tensor(&other.decompose()));
        let expected_result_2 = s.mul_tensor_identity(&other.decompose(), 37);

        assert_eq!(expected_result.size().0, 2);
        assert_eq!(expected_result.size().1, 2516);

        assert_eq!(expected_result_2.size().0, 2);
        assert_eq!(expected_result_2.size().1, 2516);

        assert_eq!(result, expected_result);
        assert_eq!(result, expected_result_2);
    }
}
