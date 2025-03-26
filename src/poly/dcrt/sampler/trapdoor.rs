#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::sync::Arc;

use crate::{
    parallel_chunk_iter, parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
        sampler::PolyTrapdoorSampler,
        Poly, PolyMatrix, PolyParams,
    },
    utils::{debug_mem, log_mem},
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{
        DCRTSquareMatTrapdoorGaussSamp, DCRTSquareMatTrapdoorGen, GetMatrixElement, Matrix,
        MatrixGen, RLWETrapdoorPair, SetMatrixElement,
    },
};

const SIGMA: f64 = 4.578;

pub struct RLWETrapdoor {
    ptr_trapdoor: Arc<UniquePtr<RLWETrapdoorPair>>,
}

pub struct DCRTMatrixPtr {
    ptr_matrix: Arc<UniquePtr<Matrix>>,
}

impl DCRTMatrixPtr {
    fn new_target_matrix(
        matrix: &DCRTPolyMatrix,
        params: &DCRTPolyParams,
        size: usize,
        target_cols: usize,
    ) -> Self {
        let mut target_matrix_ptr =
            MatrixGen(params.ring_dimension(), params.crt_depth(), params.crt_bits(), size, size);

        debug_mem(format!("target_matrix_ptr MatrixGen row={}, col={}", size, target_cols));

        for i in 0..size {
            for j in 0..target_cols {
                let entry = matrix.entry(i, j);
                let poly = entry.get_poly();
                SetMatrixElement(target_matrix_ptr.as_mut().unwrap(), i, j, poly);
            }

            if target_cols < size {
                for j in target_cols..size {
                    let zero_poly = DCRTPoly::const_zero(params);
                    let zero_poly_ptr = zero_poly.get_poly();
                    SetMatrixElement(target_matrix_ptr.as_mut().unwrap(), i, j, zero_poly_ptr);
                }
            }
        }

        Self { ptr_matrix: target_matrix_ptr.into() }
    }

    fn new_public_matrix(
        matrix: &DCRTPolyMatrix,
        n: u32,
        size: usize,
        k_res: usize,
        nrow: usize,
        ncol: usize,
    ) -> Self {
        let mut public_matrix_ptr = MatrixGen(n, size, k_res, nrow, ncol);

        debug_mem(format!("public_matrix_ptr MatrixGen row={}, col={}", nrow, ncol));

        for i in 0..nrow {
            for j in 0..ncol {
                SetMatrixElement(
                    public_matrix_ptr.as_mut().unwrap(),
                    i,
                    j,
                    matrix.entry(i, j).get_poly(),
                );
            }
        }

        Self { ptr_matrix: public_matrix_ptr.into() }
    }

    fn to_dcry_poly_matrix(
        &self,
        nrow: usize,
        ncol: usize,
        params: &DCRTPolyParams,
    ) -> DCRTPolyMatrix {
        let mut matrix_inner = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let mut row = Vec::with_capacity(ncol);
            for j in 0..ncol {
                row.push(DCRTPoly::new(GetMatrixElement(&self.ptr_matrix, i, j)));
            }
            matrix_inner.push(row);
        }

        debug_mem(format!("GetMatrixElement row={}, col={}", nrow, ncol));
        DCRTPolyMatrix::from_poly_vec(params, matrix_inner)
    }
}

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

    fn get_trapdoor_pair(&self) -> RLWETrapdoor {
        RLWETrapdoor { ptr_trapdoor: self.ptr_dcrt_trapdoor.GetTrapdoorPair().into() }
    }

    fn get_public_matrix(&self, row: usize, col: usize) -> DCRTPoly {
        DCRTPoly::new(self.ptr_dcrt_trapdoor.GetPublicMatrixElement(row, col))
    }
}

// SAFETY:
unsafe impl Send for DCRTTrapdoor {}
unsafe impl Sync for DCRTTrapdoor {}

// SAFETY:
unsafe impl Send for RLWETrapdoor {}
unsafe impl Sync for RLWETrapdoor {}

// SAFETY:
unsafe impl Send for DCRTMatrixPtr {}
unsafe impl Sync for DCRTMatrixPtr {}

pub struct DCRTPolyTrapdoorSampler {}

impl DCRTPolyTrapdoorSampler {
    pub fn new() -> Self {
        Self {}
    }
}

impl Default for DCRTPolyTrapdoorSampler {
    fn default() -> Self {
        Self::new()
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = RLWETrapdoor;
    type MatrixPtr = DCRTMatrixPtr;

    fn trapdoor(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        size: usize,
    ) -> (Self::Trapdoor, Self::M) {
        let dcrt_trapdoor = DCRTTrapdoor::new(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            SIGMA,
            2_i64,
            false,
        );
        let rlwe_trapdoor = dcrt_trapdoor.get_trapdoor_pair();
        let nrow = size;
        let ncol = (&params.modulus_bits() + 2) * size;
        let public_matrix = DCRTPolyMatrix::from_poly_vec(
            params,
            parallel_iter!(0..nrow)
                .map(|i| {
                    parallel_iter!(0..ncol).map(|j| dcrt_trapdoor.get_public_matrix(i, j)).collect()
                })
                .collect(),
        );
        (rlwe_trapdoor, public_matrix)
    }

    fn preimage(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
    ) -> Self::M {
        let size = public_matrix.row_size();
        let target_cols = target.col_size();

        assert_eq!(
            target.row_size(),
            size,
            "Target matrix should have the same number of rows as the public matrix"
        );

        // todo: real param and dummy param should have diff value
        let chunk_size = 100;
        let num_block = target_cols.div_ceil(size);
        let k = params.modulus_bits();
        debug_mem(format!(
            "preimage before loop processing with chunksize {}, out of {}",
            chunk_size, num_block
        ));

        let public_matrix = DCRTMatrixPtr::new_public_matrix(
            public_matrix,
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            (k + 2) * size,
        );

        debug_mem("SetMatrixElement public_matrix_ptr completed");
        let indices: Vec<_> = (0..num_block).collect();
        let preimages: Vec<_> = parallel_chunk_iter!(indices, chunk_size)
            .map(|chunk| {
                chunk
                    .iter()
                    .map(|&i| {
                        let start_col = i * size;
                        let end_col = (start_col + size).min(target_cols);
                        let target_block = target.slice(0, size, start_col, end_col);
                        debug_mem(format!("preimage iter: start_col = {}", start_col));
                        self.process_preimage_block(
                            params,
                            trapdoor,
                            &public_matrix,
                            &target_block,
                            size,
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .flatten()
            .collect();

        log_mem("Collected preimages");
        preimages[0].concat_columns(&preimages[1..].iter().collect::<Vec<_>>())
    }

    fn process_preimage_block(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::MatrixPtr,
        target_block: &Self::M,
        size: usize,
    ) -> Self::M {
        let n = params.ring_dimension() as usize;
        let k = params.modulus_bits();
        let target_cols = target_block.col_size();

        debug_mem(format!("Processing preimage block, target_cols={}, size={}", target_cols, size));

        let target_matrix =
            DCRTMatrixPtr::new_target_matrix(target_block, params, size, target_cols);

        debug_mem("SetMatrixElement target_matrix_ptr completed");

        let preimage_matrix = DCRTMatrixPtr {
            ptr_matrix: DCRTSquareMatTrapdoorGaussSamp(
                n as u32,
                k as u32,
                &public_matrix.ptr_matrix,
                &trapdoor.ptr_trapdoor,
                &target_matrix.ptr_matrix,
                2_i64,
                SIGMA,
            )
            .into(),
        };
        debug_mem("DCRTSquareMatTrapdoorGaussSamp completed");

        let full_preimage_matrix =
            preimage_matrix.to_dcry_poly_matrix(size * (k + 2), size, params);

        debug_mem("full_preimage_matrix generated");

        if target_cols < size {
            debug_mem("Slicing full_preimage_matrix columns");
            full_preimage_matrix.slice_columns(0, target_cols)
        } else {
            full_preimage_matrix
        }
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{sampler::DCRTPolyUniformSampler, DCRTPolyParams},
        sampler::{DistType, PolyUniformSampler},
    };

    #[test]
    fn test_trapdoor_generation() {
        let size: usize = 3;
        let sampler = DCRTPolyTrapdoorSampler::new();
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
        let size = 3;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new();
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
        let size = 4;
        let target_cols = 2;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new();
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
        let size = 4;
        let multiple = 2;
        let target_cols = size * multiple;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new();
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size and
        // target_cols is a multiple of size
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
        let size = 4;
        let target_cols = 6;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new();
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size but
        // not a multiple of size
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
