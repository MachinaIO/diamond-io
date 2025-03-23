use crate::parallel_iter;
use itertools::Itertools;
#[cfg(target_os = "linux")]
use memmap2::RemapOptions;
use memmap2::{Advice, Mmap, MmapMut, MmapOptions};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::{
    fmt::Debug,
    fs::{File, OpenOptions},
    marker::PhantomData,
    ops::{Deref, DerefMut, Index, IndexMut, Range},
    path::{Path, PathBuf},
};
use tempfile::tempfile;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Volatility {
    Persistent,
    Transient,
}

#[derive(Debug)]
pub struct MatrixStorage {
    file: File,
    nrow: usize,
    ncol: usize,
    entry_size: usize,
    volatile: Volatility,
}

impl MatrixStorage {
    pub fn new(entries: Vec<Vec<Vec<u8>>>) -> Self {
        let nrow = entries.len();
        let ncol = entries[0].len();
        debug_assert_eq!(entries.iter().all(|row| row.len() == ncol), true);
        let entry_size = entries[0][0].len();
        debug_assert_eq!(
            entries.iter().all(|row| row.iter().all(|entry| entry.len() == entry_size)),
            true
        );
        let storage = Self::new_empty(nrow, ncol, entry_size);

        let len = entry_size * nrow * ncol;
        let mut mmap = map_file_mut(&storage.file, 0, len);
        let bytes =
            entries.into_iter().flat_map(|row_vec| row_vec.into_iter().flatten()).collect_vec();
        mmap.copy_from_slice(&bytes);
        drop(mmap);
        storage
    }

    pub fn new_empty(nrow: usize, ncol: usize, entry_size: usize) -> Self {
        let len = entry_size * nrow * ncol;
        let file = tempfile().expect("failed to open file");
        file.set_len(len as u64).expect("failed to set file length");
        Self { file, nrow, ncol, entry_size, volatile: Volatility::Transient }
    }

    pub fn to_persistent(&mut self) {
        self.volatile = Volatility::Persistent;
    }

    pub fn to_transient(&mut self) {
        self.volatile = Volatility::Transient;
    }

    pub fn entries(&self, rows: Range<usize>, cols: Range<usize>) -> Vec<Vec<Vec<u8>>> {
        parallel_iter!(rows)
            .map(|i| {
                let offset = self.entry_size * (i * self.ncol + cols.start);
                let mmap = map_file(&self.file, offset, self.entry_size * cols.len());
                let row_vec = mmap.to_vec();
                let row_col_vec =
                    row_vec.chunks(self.entry_size).map(|entry| entry.to_vec()).collect();
                drop(mmap);
                row_col_vec
            })
            .collect()
    }

    pub fn replace_entries(
        &mut self,
        rows: Range<usize>,
        cols: Range<usize>,
        new_entries: Vec<Vec<Vec<u8>>>,
    ) {
        debug_assert_eq!(new_entries.len(), rows.end - rows.start);
        debug_assert_eq!(new_entries[0].len(), cols.end - cols.start);
        debug_assert_eq!(new_entries[0][0].len(), self.entry_size);

        let row_start = rows.start;
        let col_start = cols.start;
        parallel_iter!(rows).for_each(|i| {
            let offset = self.entry_size * (i * self.ncol + col_start);
            let mut mmap = map_file_mut(&self.file, offset, self.entry_size * cols.len());
            let bytes = new_entries[i - row_start].iter().cloned().flatten().collect::<Vec<_>>();
            mmap.copy_from_slice(&bytes);
            drop(mmap);
        });
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    pub fn concat_columns(&self, other: &Self) -> Self {
        debug_assert_eq!(self.nrow, other.nrow);
        debug_assert_eq!(self.entry_size, other.entry_size);
        let updated_ncol = self.ncol + other.ncol;

        let new_storage = Self::new_empty(self.nrow, updated_ncol, self.entry_size);
        parallel_iter!(0..self.nrow).for_each(|i| {
            let self_mmap = map_file(
                &self.file,
                self.entry_size * (i * self.ncol),
                self.entry_size * self.ncol,
            );
            let other_mmap = map_file(
                &other.file,
                other.entry_size * (i * other.ncol),
                other.entry_size * other.ncol,
            );
            let mut new_mmap = map_file_mut(
                &new_storage.file,
                self.entry_size * (i * updated_ncol),
                self.entry_size * updated_ncol,
            );
            new_mmap[0..self.entry_size * self.ncol].copy_from_slice(&self_mmap);
            drop(self_mmap);
            new_mmap[self.entry_size * self.ncol..].copy_from_slice(&other_mmap);
            drop(other_mmap);
            drop(new_mmap);
        });

        new_storage
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    pub fn concat_rows(&self, other: &Self) -> Self {
        debug_assert_eq!(self.ncol, other.ncol);
        debug_assert_eq!(self.entry_size, other.entry_size);
        let updated_nrow = self.nrow + other.nrow;

        let new_storage = Self::new_empty(updated_nrow, self.ncol, self.entry_size);
        parallel_iter!(0..self.nrow).for_each(|i| {
            let offset = self.entry_size * i * self.ncol;
            let self_mmap = map_file(&self.file, offset, self.entry_size * self.ncol);
            let mut new_mmap = map_file_mut(&new_storage.file, offset, self.entry_size * self.ncol);
            new_mmap.copy_from_slice(&self_mmap);
            drop(self_mmap);
            drop(new_mmap);
        });
        parallel_iter!(0..other.nrow).for_each(|i| {
            let other_mmap = map_file(
                &other.file,
                other.entry_size * i * other.ncol,
                other.entry_size * other.ncol,
            );
            let mut new_mmap = map_file_mut(
                &new_storage.file,
                self.entry_size * (self.nrow + i) * self.ncol,
                self.entry_size * self.ncol,
            );
            new_mmap.copy_from_slice(&other_mmap);
            drop(other_mmap);
            drop(new_mmap);
        });
        new_storage
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    pub fn concat_diag(&self, other: &Self) -> Self {
        debug_assert_eq!(self.entry_size, other.entry_size);
        let updated_nrow = self.nrow + other.nrow;
        let updated_ncol = self.ncol + other.ncol;

        let new_storage = Self::new_empty(updated_nrow, updated_ncol, self.entry_size);

        // Copy the first matrix (self) to the top-left corner
        parallel_iter!(0..self.nrow).for_each(|i| {
            let self_offset = self.entry_size * i * self.ncol;
            let new_offset = self.entry_size * i * updated_ncol;
            let self_mmap = map_file(&self.file, self_offset, self.entry_size * self.ncol);
            let mut new_mmap =
                map_file_mut(&new_storage.file, new_offset, self.entry_size * self.ncol);
            new_mmap.copy_from_slice(&self_mmap);
            drop(self_mmap);
            drop(new_mmap);
        });

        // Copy the second matrix (other) to the bottom-right corner
        parallel_iter!(0..other.nrow).for_each(|i| {
            let other_offset = other.entry_size * i * other.ncol;
            let new_offset = self.entry_size * ((self.nrow + i) * updated_ncol + self.ncol);
            let other_mmap = map_file(&other.file, other_offset, other.entry_size * other.ncol);
            let mut new_mmap =
                map_file_mut(&new_storage.file, new_offset, other.entry_size * other.ncol);
            new_mmap.copy_from_slice(&other_mmap);
            drop(other_mmap);
            drop(new_mmap);
        });

        new_storage
    }
}

impl Drop for MatrixStorage {
    fn drop(&mut self) {
        if self.volatile == Volatility::Transient {
            self.file.set_len(0).expect("failed to truncate file");
        }
    }
}

impl Clone for MatrixStorage {
    fn clone(&self) -> Self {
        let new_storage = Self::new_empty(self.nrow, self.ncol, self.entry_size);
        let len = self.entry_size * self.nrow * self.ncol;
        let self_mmap = map_file(&self.file, 0, len);
        let mut new_mmap = map_file_mut(&new_storage.file, 0, len);
        new_mmap.copy_from_slice(&self_mmap);
        drop(self_mmap);
        drop(new_mmap);
        new_storage
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

fn map_file_mut(file: &File, offset: usize, len: usize) -> MmapMut {
    unsafe {
        MmapOptions::new().offset(offset as u64).len(len).map_mut(file).expect("failed to map file")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a test matrix with specified dimensions
    fn create_test_matrix(nrow: usize, ncol: usize) -> MatrixStorage {
        let mut entries = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                // Create a unique value for each cell based on its position
                let value = ((i * ncol + j) as u8) % 255;
                row.push(vec![value, value + 1, value + 2, value + 3]);
            }
            entries.push(row);
        }
        MatrixStorage::new(entries)
    }

    #[test]
    fn test_matrix_storage_new_and_entries() {
        // Create a 3x4 matrix
        let nrow = 3;
        let ncol = 4;
        let entry_size = 4; // Each entry is 4 bytes

        // Create test data
        let mut entries = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                let value = ((i * ncol + j) as u8) % 255;
                row.push(vec![value, value + 1, value + 2, value + 3]);
            }
            entries.push(row);
        }

        // Create matrix
        let matrix = MatrixStorage::new(entries.clone());

        // Verify dimensions
        assert_eq!(matrix.nrow, nrow);
        assert_eq!(matrix.ncol, ncol);
        assert_eq!(matrix.entry_size, entry_size);

        // Retrieve all entries and verify they match the original data
        let retrieved_entries = matrix.entries(0..nrow, 0..ncol);
        assert_eq!(retrieved_entries, entries);
    }

    #[test]
    fn test_matrix_storage_replace_entries() {
        // Create a 3x4 matrix
        let matrix = create_test_matrix(3, 4);

        // Create new entries for a 2x2 submatrix
        let new_entries = vec![
            vec![vec![100, 101, 102, 103], vec![110, 111, 112, 113]],
            vec![vec![120, 121, 122, 123], vec![130, 131, 132, 133]],
        ];

        // Replace a 2x2 submatrix starting at position (1, 1)
        let mut modified_matrix = matrix.clone();
        modified_matrix.replace_entries(1..3, 1..3, new_entries.clone());

        // Retrieve the modified submatrix and verify it matches the new entries
        let retrieved_entries = modified_matrix.entries(1..3, 1..3);
        assert_eq!(retrieved_entries, new_entries);

        // Verify that the rest of the matrix is unchanged
        let top_left = modified_matrix.entries(0..1, 0..1);
        assert_eq!(top_left, matrix.entries(0..1, 0..1));

        let top_right = modified_matrix.entries(0..1, 3..4);
        assert_eq!(top_right, matrix.entries(0..1, 3..4));

        let bottom_left = modified_matrix.entries(2..3, 0..1);
        assert_eq!(bottom_left, matrix.entries(2..3, 0..1));
    }

    #[test]
    fn test_matrix_storage_concat_columns() {
        // Create two matrices with the same number of rows
        let matrix1 = create_test_matrix(3, 2);
        let matrix2 = create_test_matrix(3, 3);

        // Concatenate the matrices horizontally
        let result = matrix1.concat_columns(&matrix2);

        // Verify dimensions
        assert_eq!(result.nrow, 3);
        assert_eq!(result.ncol, 5);
        assert_eq!(result.entry_size, 4);

        // Verify the left part matches matrix1
        let left_part = result.entries(0..3, 0..2);
        assert_eq!(left_part, matrix1.entries(0..3, 0..2));

        // Verify the right part matches matrix2
        let right_part = result.entries(0..3, 2..5);
        assert_eq!(right_part, matrix2.entries(0..3, 0..3));
    }

    #[test]
    fn test_matrix_storage_concat_rows() {
        // Create two matrices with the same number of columns
        let matrix1 = create_test_matrix(2, 3);
        let matrix2 = create_test_matrix(3, 3);

        // Concatenate the matrices vertically
        let result = matrix1.concat_rows(&matrix2);

        // Verify dimensions
        assert_eq!(result.nrow, 5);
        assert_eq!(result.ncol, 3);
        assert_eq!(result.entry_size, 4);

        // Verify the top part matches matrix1
        let top_part = result.entries(0..2, 0..3);
        assert_eq!(top_part, matrix1.entries(0..2, 0..3));

        // Verify the bottom part matches matrix2
        let bottom_part = result.entries(2..5, 0..3);
        assert_eq!(bottom_part, matrix2.entries(0..3, 0..3));
    }

    #[test]
    fn test_matrix_storage_concat_diag() {
        // Create two matrices
        let matrix1 = create_test_matrix(2, 3);
        let matrix2 = create_test_matrix(3, 2);

        // Concatenate the matrices diagonally
        let result = matrix1.concat_diag(&matrix2);

        // Verify dimensions
        assert_eq!(result.nrow, 5);
        assert_eq!(result.ncol, 5);
        assert_eq!(result.entry_size, 4);

        // Verify the top-left part matches matrix1
        let top_left = result.entries(0..2, 0..3);
        assert_eq!(top_left, matrix1.entries(0..2, 0..3));

        // Verify the bottom-right part matches matrix2
        let bottom_right = result.entries(2..5, 3..5);
        assert_eq!(bottom_right, matrix2.entries(0..3, 0..2));

        // Verify the top-right part is zeros (or empty)
        let top_right = result.entries(0..2, 3..5);
        for row in top_right {
            for entry in row {
                assert_eq!(entry, vec![0, 0, 0, 0]);
            }
        }

        // Verify the bottom-left part is zeros (or empty)
        let bottom_left = result.entries(2..5, 0..3);
        for row in bottom_left {
            for entry in row {
                assert_eq!(entry, vec![0, 0, 0, 0]);
            }
        }
    }

    #[test]
    fn test_matrix_storage_volatility() {
        // Create a matrix with default transient volatility
        let mut matrix = create_test_matrix(2, 2);
        assert_eq!(matrix.volatile, Volatility::Transient);

        // Change to persistent
        matrix.to_persistent();
        assert_eq!(matrix.volatile, Volatility::Persistent);

        // Change back to transient
        matrix.to_transient();
        assert_eq!(matrix.volatile, Volatility::Transient);
    }
}
