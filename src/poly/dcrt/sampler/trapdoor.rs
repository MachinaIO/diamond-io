use std::sync::Arc;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
    sampler::PolyTrapdoorSampler,
    PolyMatrix, PolyParams,
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{
        DCRTSquareMatGaussSamp, DCRTSquareMatTrapdoorGen, GetMatrixElement, MatrixGen,
        RLWETrapdoorPair, SetMatrixElement,
    },
};

pub struct DCRTPolyTrapdoorSampler {
    params: DCRTPolyParams,
    base: usize,
    sigma: f64,
    d: usize,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(params: DCRTPolyParams, base: usize, sigma: f64, d: usize) -> Self {
        Self { params, base, sigma, d }
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = Arc<UniquePtr<RLWETrapdoorPair>>;

    fn trapdoor(&self) -> (Self::Trapdoor, Self::M) {
        let trapdoor_output = DCRTSquareMatTrapdoorGen(
            self.params.ring_dimension(),
            self.params.size(),
            self.params.k_res(),
            self.d,
            self.sigma,
            self.base as i64,
            false,
        );
        let trapdoor = trapdoor_output.GetTrapdoorPair();
        let nrow = self.d;
        let ncol = (&self.params.modulus_bits() + 2) * self.d;

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

        let public_matrix = DCRTPolyMatrix::from_poly_vec(&self.params, matrix_inner);

        (trapdoor.into(), public_matrix)
    }

    fn preimage(
        &self,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
    ) -> Self::M {
        let n = self.params.ring_dimension() as usize;
        let k = self.params.modulus_bits();

        let mut public_matrix_ptr = MatrixGen(
            self.params.ring_dimension(),
            self.params.size(),
            self.params.k_res(),
            self.d,
            (k + 2) * self.d,
        );

        for i in 0..self.d {
            for j in 0..(k + 2) * self.d {
                let poly = public_matrix.entry(i, j).get_poly();
                SetMatrixElement(public_matrix_ptr.as_mut().unwrap(), i, j, poly);
            }
        }

        let mut target_matrix_ptr = MatrixGen(
            self.params.ring_dimension(),
            self.params.size(),
            self.params.k_res(),
            self.d,
            self.d,
        );

        for i in 0..self.d {
            for j in 0..self.d {
                let poly = target.entry(i, j).get_poly();
                SetMatrixElement(target_matrix_ptr.as_mut().unwrap(), i, j, poly);
            }
        }

        let preimage_matrix_ptr = DCRTSquareMatGaussSamp(
            n as u32,
            k as u32,
            &public_matrix_ptr,
            trapdoor,
            &target_matrix_ptr,
            self.base as i64,
            self.sigma,
        );

        let nrow = self.d * (k + 2);
        let ncol = self.d;

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

        DCRTPolyMatrix::from_poly_vec(&self.params, matrix_inner)
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
        let params = DCRTPolyParams::default();
        let base = 2;
        let sigma = 4.57825;
        let d = 3;
        let sampler = DCRTPolyTrapdoorSampler::new(params, base, sigma, d);

        let (_, public_matrix) = sampler.trapdoor();

        let expected_rows = d;
        let expected_cols = (&sampler.params.modulus_bits() + 2) * d;

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
        let d = 3;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(params.clone(), base, sigma, d);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor();

        let uniform_sampler = DCRTPolyUniformSampler::new(params.clone());
        let target = uniform_sampler.sample_uniform(d, d, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&trapdoor, &public_matrix, &target);

        let expected_rows = d * (k + 2);
        let expected_cols = d;

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
