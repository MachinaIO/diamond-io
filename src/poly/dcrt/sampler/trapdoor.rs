#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::sync::Arc;

use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix},
        sampler::PolyTrapdoorSampler,
        Poly, PolyMatrix, PolyParams,
    },
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{
        DCRTSquareMatTrapdoorGaussSamp, DCRTSquareMatTrapdoorGen, GetMatrixElement, MatrixGen,
        SetMatrixElement,
    },
};

pub struct DCRTTrapdoor {
    ptr_dcrt_trapdoor: Arc<UniquePtr<openfhe::ffi::DCRTTrapdoor>>,
}

impl DCRTTrapdoor {
    fn new(
        n: u32,
        size: usize,
        k_res: usize,
        d: usize,
        sigma: f64,
        base: i64,
        balanced: bool,
    ) -> Self {
        let ptr_dcrt_trapdoor = DCRTSquareMatTrapdoorGen(n, size, k_res, d, sigma, base, balanced);
        Self { ptr_dcrt_trapdoor: ptr_dcrt_trapdoor.into() }
    }

    fn get_public_matrix_element(&self, row: usize, col: usize) -> DCRTPoly {
        DCRTPoly::new(GetMatrixElement(&self.ptr_dcrt_trapdoor.GetPublicMatrix(), row, col))
    }

    fn get_trapdoor_first_element(&self, row: usize, col: usize) -> DCRTPoly {
        DCRTPoly::new(GetMatrixElement(&self.ptr_dcrt_trapdoor.GetTrapdoorFirst(), row, col))
    }

    fn get_trapdoor_second_element(&self, row: usize, col: usize) -> DCRTPoly {
        DCRTPoly::new(GetMatrixElement(&self.ptr_dcrt_trapdoor.GetTrapdoorSecond(), row, col))
    }
}

// SAFETY:
unsafe impl Send for DCRTTrapdoor {}
unsafe impl Sync for DCRTTrapdoor {}

pub struct DCRTPolyTrapdoorSampler {
    base: usize,
    sigma: f64,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(base: usize, sigma: f64) -> Self {
        Self { base, sigma }
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;

    fn trapdoor(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        size: usize,
    ) -> (Self::M, Self::M, Self::M) {
        let dcrt_trapdoor = DCRTTrapdoor::new(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            self.sigma,
            self.base as i64,
            false,
        );

        let nrow_public_mat = size;
        let ncol_public_mat = (&params.modulus_bits() + 2) * size;
        let public_matrix = DCRTPolyMatrix::from_poly_vec(
            params,
            parallel_iter!(0..nrow_public_mat)
                .map(|i| {
                    parallel_iter!(0..ncol_public_mat)
                        .map(|j| dcrt_trapdoor.get_public_matrix_element(i, j))
                        .collect()
                })
                .collect(),
        );

        let nrow_trapdoor = size;
        let ncol_trapdoor = params.modulus_bits() * size;
        let trapdoor_first = DCRTPolyMatrix::from_poly_vec(
            params,
            parallel_iter!(0..nrow_trapdoor)
                .map(|i| {
                    parallel_iter!(0..ncol_trapdoor)
                        .map(|j| dcrt_trapdoor.get_trapdoor_first_element(i, j))
                        .collect()
                })
                .collect(),
        );
        let trapdoor_second = DCRTPolyMatrix::from_poly_vec(
            params,
            parallel_iter!(0..nrow_trapdoor)
                .map(|i| {
                    parallel_iter!(0..ncol_trapdoor)
                        .map(|j| dcrt_trapdoor.get_trapdoor_second_element(i, j))
                        .collect()
                })
                .collect(),
        );

        (public_matrix, trapdoor_first, trapdoor_second)
    }

    fn preimage(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        public_matrix: &Self::M,
        trapdoor_first: &Self::M,
        trapdoor_second: &Self::M,
        target: &Self::M,
    ) -> Self::M {
        let n = params.ring_dimension() as usize;
        let k = params.modulus_bits();
        let size = public_matrix.row_size();
        let target_cols = target.col_size();

        assert_eq!(
            target.row_size(),
            size,
            "Target matrix should have the same number of rows as the public matrix"
        );

        // Case 1: Target columns is greater than size
        if target_cols > size {
            let full_blocks = target_cols / size;
            let remaining_cols = target_cols % size;
            let total_blocks = if remaining_cols > 0 { full_blocks + 1 } else { full_blocks };

            let preimages: Vec<_> = parallel_iter!(0..total_blocks)
                .map(|block| {
                    let start_col = block * size;

                    // Calculate end_col based on whether this is the last block with remaining
                    // columns
                    let end_col = if block == full_blocks && remaining_cols > 0 {
                        start_col + remaining_cols
                    } else {
                        start_col + size
                    };

                    // Process the block
                    let target_block = target.slice(0, size, start_col, end_col);
                    self.preimage(
                        params,
                        public_matrix,
                        trapdoor_first,
                        trapdoor_second,
                        &target_block,
                    )
                })
                .collect();

            // Concatenate all preimages horizontally
            return preimages[0].concat_columns(&preimages[1..].iter().collect::<Vec<_>>());
        }

        // Case 2: Target columns is equal or less than size
        let mut public_matrix_ptr = MatrixGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            (k + 2) * size,
        );

        for i in 0..size {
            for j in 0..(k + 2) * size {
                let poly = public_matrix.entry(i, j).get_poly();
                SetMatrixElement(public_matrix_ptr.as_mut().unwrap(), i, j, poly);
            }
        }

        let mut trapdoor_first_ptr = MatrixGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            k * size,
        );

        let mut trapdoor_second_ptr = MatrixGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            k * size,
        );

        for i in 0..size {
            for j in 0..k * size {
                let poly_first = trapdoor_first.entry(i, j).get_poly();
                let poly_second = trapdoor_second.entry(i, j).get_poly();
                SetMatrixElement(trapdoor_first_ptr.as_mut().unwrap(), i, j, poly_first);
                SetMatrixElement(trapdoor_second_ptr.as_mut().unwrap(), i, j, poly_second);
            }
        }

        let mut target_matrix_ptr =
            MatrixGen(params.ring_dimension(), params.crt_depth(), params.crt_bits(), size, size);

        for i in 0..size {
            for j in 0..target_cols {
                let poly = target.entry(i, j).get_poly();
                SetMatrixElement(target_matrix_ptr.as_mut().unwrap(), i, j, poly);
            }

            // Pad the remaining columns with zeros if target_cols < size
            if target_cols < size {
                for j in target_cols..size {
                    let zero_poly = DCRTPoly::const_zero(params);
                    let zero_poly_ptr = zero_poly.get_poly();
                    SetMatrixElement(target_matrix_ptr.as_mut().unwrap(), i, j, zero_poly_ptr);
                }
            }
        }

        let preimage_matrix_ptr = DCRTSquareMatTrapdoorGaussSamp(
            n as u32,
            k as u32,
            &public_matrix_ptr,
            &trapdoor_first_ptr,
            &trapdoor_second_ptr,
            &target_matrix_ptr,
            self.base as i64,
            self.sigma,
        );

        let nrow = size * (k + 2);
        let ncol = size;

        let mut matrix_inner = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                let poly = GetMatrixElement(&preimage_matrix_ptr, i, j);
                let dcrt_poly = DCRTPoly::new(poly);
                row.push(dcrt_poly);
            }
            matrix_inner.push(row);
        }

        let full_preimage = DCRTPolyMatrix::from_poly_vec(params, matrix_inner);

        // If the target matrix has fewer columns than size, slice the preimage matrix
        if target_cols < size {
            full_preimage.slice_columns(0, target_cols)
        } else {
            full_preimage
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{sampler::DCRTPolyUniformSampler, DCRTPolyParams},
        sampler::{DistType, PolyUniformSampler},
    };

    #[test]
    fn test_trapdoor_generation() {
        let base = 2;
        let sigma = 4.57825;
        let size: usize = 3;
        let sampler = DCRTPolyTrapdoorSampler::new(base, sigma);
        let params = DCRTPolyParams::default();

        let (public_matrix, _, _) = sampler.trapdoor(&params, size);

        let expected_rows = size;
        let expected_cols = (&params.modulus_bits() + 2) * size;

        assert_eq!(
            public_matrix.row_size(),
            expected_rows,
            "Public matrix should have the correct number of rows"
        );
        assert_eq!(
            public_matrix.col_size(),
            expected_cols,
            "Public matrix should have the correct number of columns"
        );

        // Verify that all entries in the matrix are valid DCRTPolys
        for i in 0..public_matrix.row_size() {
            for j in 0..public_matrix.col_size() {
                let poly = public_matrix.entry(i, j);
                assert!(!poly.get_poly().is_null(), "Matrix entry should be a valid DCRTPoly");
            }
        }
    }

    #[test]
    fn test_preimage_generation() {
        let params = DCRTPolyParams::default();
        let base = 2;
        let sigma = 4.57825;
        let size = 3;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(base, sigma);
        let (public_matrix, trapdoor_first, trapdoor_second) =
            trapdoor_sampler.trapdoor(&params, size);

        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, size, size, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(
            &params,
            &public_matrix,
            &trapdoor_first,
            &trapdoor_second,
            &target,
        );

        let expected_rows = size * (k + 2);
        let expected_cols = size;

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns"
        );

        // public_matrix * preimage should be equal to target
        let product = public_matrix * &preimage;
        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }

    #[test]
    fn test_preimage_generation_non_square_target_lt() {
        let params = DCRTPolyParams::default();
        let base = 2;
        let sigma = 4.57825;
        let size = 4;
        let target_cols = 2;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(base, sigma);
        let (public_matrix, trapdoor_first, trapdoor_second) =
            trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols < size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        // Compute the preimage
        let preimage = trapdoor_sampler.preimage(
            &params,
            &public_matrix,
            &trapdoor_first,
            &trapdoor_second,
            &target,
        );

        // Verify dimensions of the preimage matrix
        let expected_rows = size * (k + 2);
        let expected_cols = target_cols; // Preimage should be sliced to match target columns

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns (sliced to match target)"
        );

        // Verify that public_matrix * preimage = target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }

    #[test]
    fn test_preimage_generation_non_square_target_gt_multiple() {
        let params = DCRTPolyParams::default();
        let base = 2;
        let sigma = 4.57825;
        let size = 4;
        let multiple = 2;
        let target_cols = size * multiple;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(base, sigma);
        let (public_matrix, trapdoor_first, trapdoor_second) =
            trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size and
        // target_cols is a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        // Compute the preimage
        let preimage = trapdoor_sampler.preimage(
            &params,
            &public_matrix,
            &trapdoor_first,
            &trapdoor_second,
            &target,
        );

        // Verify dimensions of the preimage matrix
        let expected_rows = size * (k + 2);
        let expected_cols = target_cols;

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns (equal to target columns)"
        );

        // Verify that public_matrix * preimage = target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }

    #[test]
    fn test_preimage_generation_non_square_target_gt_non_multiple() {
        let params = DCRTPolyParams::default();
        let base = 2;
        let sigma = 4.57825;
        let size = 4;
        let target_cols = 6;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(base, sigma);
        let (public_matrix, trapdoor_first, trapdoor_second) =
            trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size but
        // not a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        // Compute the preimage
        let preimage = trapdoor_sampler.preimage(
            &params,
            &public_matrix,
            &trapdoor_first,
            &trapdoor_second,
            &target,
        );

        // Verify dimensions of the preimage matrix
        let expected_rows = size * (k + 2);
        let expected_cols = target_cols;

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns (equal to target columns)"
        );

        // Verify that public_matrix * preimage = target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }
}
