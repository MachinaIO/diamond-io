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

    fn get_trapdoor_pair(&self) -> UniquePtr<RLWETrapdoorPair> {
        self.ptr_dcrt_trapdoor.GetTrapdoorPair()
    }

    fn get_public_matrix(&self) -> UniquePtr<Matrix> {
        self.ptr_dcrt_trapdoor.GetPublicMatrix()
    }
}

// SAFETY:
unsafe impl Send for DCRTTrapdoor {}
unsafe impl Sync for DCRTTrapdoor {}

pub struct DCRTPolyTrapdoorSampler {
    pubic_matrix_ptr: UniquePtr<Matrix>,
    trapdoor_pair_ptr: UniquePtr<RLWETrapdoorPair>,
    size: usize,
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;

    fn new(params: &<<Self::M as PolyMatrix>::P as Poly>::Params, size: usize) -> Self {
        let dcrt_trapdoor = DCRTTrapdoor::new(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            size,
            SIGMA,
            2_i64,
            false,
        );
        let trapdoor_pair_ptr = dcrt_trapdoor.get_trapdoor_pair();
        let pubic_matrix_ptr = dcrt_trapdoor.get_public_matrix();

        Self { pubic_matrix_ptr, trapdoor_pair_ptr, size }
    }

    fn preimage(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        target: &Self::M,
    ) -> UniquePtr<Matrix> {
        let target_cols = target.col_size();
        let size = self.size;
        assert_eq!(
            target.row_size(),
            size,
            "Target matrix should have the same number of rows as the public matrix"
        );

        debug_mem("preimage before loop processing");

        let num_block = target_cols.div_ceil(size);
        let preimages: Vec<UniquePtr<Matrix>> = parallel_iter!(0..num_block)
            .map(|i| {
                let start_col = i * size;
                let end_col = (start_col + size).min(target_cols);
                let target_block = target.slice(0, size, start_col, end_col);
                debug_mem(format!("preimage iter : start_col = {}", start_col));

                self.process_preimage_block(params, &target_block)
            })
            .collect();

        log_mem("Collected preimages");
        // preimages[0].concat_columns(&preimages[1..].iter().collect::<Vec<_>>())
    }

    fn process_preimage_block(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        target_block: &Self::M,
    ) -> UniquePtr<Matrix> {
        let size = self.size;
        let n = params.ring_dimension() as usize;
        let k = params.modulus_bits();
        // let size = public_matrix.row_size();
        let target_cols = target_block.col_size();

        debug_mem("Processing preimage block");

        debug_mem("SetMatrixElement public_matrix_ptr completed");

        let mut target_matrix_ptr =
            MatrixGen(params.ring_dimension(), params.crt_depth(), params.crt_bits(), size, size);

        debug_mem("target_matrix_ptr generated");

        for i in 0..size {
            for j in 0..target_cols {
                let entry = target_block.entry(i, j);
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

        debug_mem("SetMatrixElement target_matrix_ptr completed");

        let preimage_matrix_ptr = DCRTSquareMatTrapdoorGaussSamp(
            n as u32,
            k as u32,
            &self.pubic_matrix_ptr,
            &self.trapdoor_pair_ptr,
            &target_matrix_ptr,
            2_i64,
            SIGMA,
        );

        preimage_matrix_ptr
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
        let params = DCRTPolyParams::default();

        let sampler = DCRTPolyTrapdoorSampler::new(&params, size);

        let expected_rows = size;
        let expected_cols = (&params.modulus_bits() + 2) * size;

        // assert_eq!(
        //     public_matrix.row_size(),
        //     expected_rows,
        //     "Public matrix should have the correct number of rows"
        // );
        // assert_eq!(
        //     public_matrix.col_size(),
        //     expected_cols,
        //     "Public matrix should have the correct number of columns"
        // );

        // // Verify that all entries in the matrix are valid DCRTPolys
        // for i in 0..public_matrix.row_size() {
        //     for j in 0..public_matrix.col_size() {
        //         let poly = public_matrix.entry(i, j);
        //         assert!(!poly.get_poly().is_null(), "Matrix entry should be a valid DCRTPoly");
        //     }
        // }
    }

    #[test]
    fn test_preimage_generation() {
        let params = DCRTPolyParams::default();
        let size = 3;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, size);

        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, size, size, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&params, &target);

        let expected_rows = size * (k + 2);
        let expected_cols = size;

        // assert_eq!(
        //     preimage.row_size(),
        //     expected_rows,
        //     "Preimage matrix should have the correct number of rows"
        // );

        // assert_eq!(
        //     preimage.col_size(),
        //     expected_cols,
        //     "Preimage matrix should have the correct number of columns"
        // );

        // public_matrix * preimage should be equal to target
        // todo
        // let product = trapdoor_sampler.pubic_matrix_ptr * &preimage;
        // assert_eq!(product, target, "Product of public matrix and preimage should equal target");
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
