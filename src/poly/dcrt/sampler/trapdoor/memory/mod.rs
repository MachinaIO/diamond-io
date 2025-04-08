use rayon::iter::ParallelIterator;
use std::sync::Arc;

use crate::{
    parallel_iter,
    poly::{
        dcrt::{cpp_matrix::CppMatrix, DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
        sampler::PolyTrapdoorSampler,
        Poly, PolyMatrix, PolyParams,
    },
    utils::log_mem,
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{DCRTSquareMatTrapdoorGaussSamp, DCRTSquareMatTrapdoorGen, RLWETrapdoorPair},
};

pub struct RLWETrapdoor {
    ptr_trapdoor: Arc<UniquePtr<RLWETrapdoorPair>>,
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

    fn get_public_matrix(&self) -> CppMatrix {
        CppMatrix::new(self.ptr_dcrt_trapdoor.GetPublicMatrix())
    }
}

// SAFETY:
unsafe impl Send for DCRTTrapdoor {}
unsafe impl Sync for DCRTTrapdoor {}

// SAFETY:
unsafe impl Send for RLWETrapdoor {}
unsafe impl Sync for RLWETrapdoor {}

pub struct DCRTPolyTrapdoorSampler {
    sigma: f64,
    base: u32,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(params: &DCRTPolyParams, sigma: f64) -> Self {
        let base = 1 << params.base_bits();
        Self { sigma, base }
    }
}

impl DCRTPolyTrapdoorSampler {
    fn process_preimage_block(
        &self,
        params: &DCRTPolyParams,
        trapdoor: &RLWETrapdoor,
        public_matrix: &CppMatrix,
        target_block: &DCRTPolyMatrix,
        size: usize,
    ) -> DCRTPolyMatrix {
        let n = params.ring_dimension() as usize;
        let k = params.modulus_bits();
        let target_cols = target_block.col_size();
        log_mem(format!("Processing preimage block, target_cols={}, size={}", target_cols, size));
        let target_matrix = target_block.to_cpp_matrix_ptr();
        log_mem("SetMatrixElement target_matrix_ptr completed");

        let preimage_matrix = CppMatrix::new(DCRTSquareMatTrapdoorGaussSamp(
            n as u32,
            k as u32,
            &public_matrix.inner,
            &trapdoor.ptr_trapdoor,
            &target_matrix.inner,
            self.base.into(),
            self.sigma,
        ));
        log_mem("DCRTSquareMatTrapdoorGaussSamp completed");

        let full_preimage_matrix = DCRTPolyMatrix::from_cpp_matrix_ptr(params, &preimage_matrix);
        log_mem("full_preimage_matrix generated");

        if target_cols < size {
            log_mem("Slicing full_preimage_matrix columns");
            full_preimage_matrix.slice_columns(0, target_cols)
        } else {
            full_preimage_matrix
        }
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = RLWETrapdoor;

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
            self.sigma,
            self.base.into(),
            false,
        );
        let rlwe_trapdoor = dcrt_trapdoor.get_trapdoor_pair();
        let nrow = size;
        let ncol = (&params.modulus_digits() + 2) * size;
        let public_matrix =
            DCRTPolyMatrix::from_cpp_matrix_ptr(&params, &dcrt_trapdoor.get_public_matrix());
        log_mem(format!(
            "public matrix: {} {} || {} {}",
            public_matrix.col_size(),
            public_matrix.row_size(),
            nrow,
            ncol
        ));
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

        let num_block = target_cols.div_ceil(size);
        log_mem(format!("preimage before loop processing out of {}", num_block));
        let public_matrix = public_matrix.to_cpp_matrix_ptr();

        log_mem("SetMatrixElement public_matrix_ptr completed");

        let preimages: Vec<_> = parallel_iter!(0..num_block)
            .map(|i| {
                let start_col = i * size;
                let end_col = (start_col + size).min(target_cols);
                let target_block = target.slice(0, size, start_col, end_col);
                log_mem(format!("preimage iter : start_col = {}", start_col));

                self.process_preimage_block(params, trapdoor, &public_matrix, &target_block, size)
            })
            .collect();

        log_mem("Collected preimages");
        preimages[0].concat_columns(&preimages[1..].iter().collect::<Vec<_>>())
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::{
        poly::{
            dcrt::{sampler::DCRTPolyUniformSampler, DCRTPolyParams},
            sampler::{DistType, PolyUniformSampler},
        },
        utils::init_tracing,
    };

    const SIGMA: f64 = 4.578;

    #[test]
    fn test_trapdoor_generation() {
        let size: usize = 3;

        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);

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
        let size = 4;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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

    #[test]
    fn test_preimage_generation_base_8() {
        init_tracing();
        let params = DCRTPolyParams::new(4, 2, 17, 3);
        let size = 4;
        let target_cols = 6;
        let k = params.modulus_digits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        log_mem(format!(
            "public_matrix :{} {}",
            public_matrix.col_size(),
            public_matrix.row_size()
        ));

        // Create a non-square target matrix (size x target_cols) such that target_cols > size
        // target_cols is not a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);
        log_mem(format!("target :{} {}", target.col_size(), target.row_size()));

        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);
        log_mem(format!("preimage :{} {}", preimage.col_size(), preimage.row_size()));
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

        // public_matrix * preimage should be equal to target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }
}
