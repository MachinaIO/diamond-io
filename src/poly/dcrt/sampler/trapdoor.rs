use std::sync::Arc;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix},
    sampler::PolyTrapdoorSampler,
    Poly, PolyMatrix, PolyParams,
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{
        DCRTSquareMatTrapdoorGaussSamp, DCRTSquareMatTrapdoorGen, GetMatrixElement, MatrixGen,
        RLWETrapdoorPair, SetMatrixElement,
    },
};

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
    type Trapdoor = Arc<UniquePtr<RLWETrapdoorPair>>;

    fn trapdoor(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        size: usize,
    ) -> (Self::Trapdoor, Self::M) {
        let trapdoor_output = DCRTSquareMatTrapdoorGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            self.sigma,
            self.base as i64,
            false,
        );
        let trapdoor = trapdoor_output.GetTrapdoorPair();
        let nrow = size;
        let ncol = (&params.modulus_bits() + 2) * size;

        // Construct the public matrix from its elements
        let mut matrix_inner = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                let poly = trapdoor_output.GetPublicMatrixElement(i, j);
                let dcrt_poly = DCRTPoly::new(poly);
                row.push(dcrt_poly);
            }
            matrix_inner.push(row);
        }

        let public_matrix = DCRTPolyMatrix::from_poly_vec(params, matrix_inner);

        (trapdoor.into(), public_matrix)
    }

    fn preimage(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
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
            // Calculate how many full blocks of size columns we have
            let full_blocks = target_cols / size;
            // Calculate the remaining columns (could be 0)
            let remaining_cols = target_cols % size;

            // Total number of blocks (full blocks + 1 partial block if there are remaining columns)
            let total_blocks = if remaining_cols > 0 { full_blocks + 1 } else { full_blocks };

            let mut preimages = Vec::with_capacity(total_blocks);

            // Process each block
            for block in 0..total_blocks {
                let start_col = block * size;

                // For the last block with remaining columns
                if block == full_blocks && remaining_cols > 0 {
                    let end_col = start_col + remaining_cols;
                    let target_block = target.slice(0, size, start_col, end_col);
                    let preimage_block =
                        self.preimage(params, trapdoor, public_matrix, &target_block);
                    preimages.push(preimage_block);
                } else if block < full_blocks {
                    // For full blocks
                    let end_col = start_col + size;
                    let target_block = target.slice(0, size, start_col, end_col);
                    let preimage_block =
                        self.preimage(params, trapdoor, public_matrix, &target_block);
                    preimages.push(preimage_block);
                }
            }

            // Concatenate all preimages horizontally
            if preimages.len() == 1 {
                return preimages[0].clone();
            } else {
                let first = preimages[0].clone();
                let rest = &preimages[1..];
                return first.concat_columns(rest);
            }
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
            trapdoor,
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

        let (_, public_matrix) = sampler.trapdoor(&params, size);

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
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, size, size, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

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
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols < size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        // Compute the preimage
        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

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
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size and target_cols is a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        // Compute the preimage
        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

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
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size but not a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        // Compute the preimage
        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

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
