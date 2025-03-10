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
    /// size of the target square matrix (size X size)
    size: usize,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(base: usize, sigma: f64, size: usize) -> Self {
        Self { base, sigma, size }
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = Arc<UniquePtr<RLWETrapdoorPair>>;

    fn trapdoor(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
    ) -> (Self::Trapdoor, Self::M) {
        let trapdoor_output = DCRTSquareMatTrapdoorGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            self.size,
            self.sigma,
            self.base as i64,
            false,
        );
        let trapdoor = trapdoor_output.GetTrapdoorPair();
        let nrow = self.size;
        let ncol = (&params.modulus_bits() + 2) * self.size;

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

        let mut public_matrix_ptr = MatrixGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            self.size,
            (k + 2) * self.size,
        );

        for i in 0..self.size {
            for j in 0..(k + 2) * self.size {
                let poly = public_matrix.entry(i, j).get_poly();
                SetMatrixElement(public_matrix_ptr.as_mut().unwrap(), i, j, poly);
            }
        }

        let mut target_matrix_ptr = MatrixGen(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            self.size,
            self.size,
        );

        // if target.col_size() > self.size {
        //     panic!("Target matrix should have the same size as the public matrix");
        // }

        for i in 0..self.size {
            for j in 0..self.size {
                let poly = target.entry(i, j).get_poly();
                SetMatrixElement(target_matrix_ptr.as_mut().unwrap(), i, j, poly);
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

        let nrow = self.size * (k + 2);
        let ncol = self.size;

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

        DCRTPolyMatrix::from_poly_vec(params, matrix_inner)
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
        let sampler = DCRTPolyTrapdoorSampler::new(base, sigma, size);
        let params = DCRTPolyParams::default();

        let (_, public_matrix) = sampler.trapdoor(&params);

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
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(base, sigma, size);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params);

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
}
