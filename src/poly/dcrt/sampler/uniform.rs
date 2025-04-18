use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix},
        sampler::{DistType, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use openfhe::ffi;
use rayon::prelude::*;
#[cfg(feature = "disk")]
use std::ops::Range;

pub struct DCRTPolyUniformSampler {}

impl Default for DCRTPolyUniformSampler {
    fn default() -> Self {
        Self::new()
    }
}

impl PolyUniformSampler for DCRTPolyUniformSampler {
    type M = DCRTPolyMatrix;

    fn new() -> Self {
        Self {}
    }

    fn sample_poly(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        dist: &DistType,
    ) -> <Self::M as PolyMatrix>::P {
        let sampled_poly = match dist {
            DistType::FinRingDist => ffi::DCRTPolyGenFromDug(
                params.ring_dimension(),
                params.crt_depth(),
                params.crt_bits(),
            ),
            DistType::GaussDist { sigma } => ffi::DCRTPolyGenFromDgg(
                params.ring_dimension(),
                params.crt_depth(),
                params.crt_bits(),
                *sigma,
            ),
            DistType::BitDist => ffi::DCRTPolyGenFromBug(
                params.ring_dimension(),
                params.crt_depth(),
                params.crt_bits(),
            ),
        };
        if sampled_poly.is_null() {
            panic!("Attempted to dereference a null pointer");
        }
        DCRTPoly::new(sampled_poly)
    }

    fn sample_uniform(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        nrow: usize,
        ncol: usize,
        dist: DistType,
    ) -> Self::M {
        #[cfg(feature = "disk")]
        {
            let mut new_matrix = DCRTPolyMatrix::new_empty(params, nrow, ncol);
            let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
                parallel_iter!(row_offsets)
                    .map(|_| {
                        parallel_iter!(col_offsets.clone())
                            .map(|_| self.sample_poly(params, &dist))
                            .collect()
                    })
                    .collect()
            };
            new_matrix.replace_entries(0..nrow, 0..ncol, f);
            new_matrix
        }
        #[cfg(not(feature = "disk"))]
        {
            let c: Vec<Vec<DCRTPoly>> = parallel_iter!(0..nrow)
                .map(|_| {
                    parallel_iter!(0..ncol)
                        .map(|_| {
                            let sampled_poly = self.sample_poly(params, &dist);
                            if sampled_poly.get_poly().is_null() {
                                panic!("Attempted to dereference a null pointer");
                            }
                            sampled_poly
                        })
                        .collect()
                })
                .collect();

            DCRTPolyMatrix::from_poly_vec(params, c)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;

    #[test]
    fn test_ring_dist() {
        let params = DCRTPolyParams::default();

        // Test FinRingDist
        let sampler = DCRTPolyUniformSampler::new();
        let matrix1 = sampler.sample_uniform(&params, 20, 5, DistType::FinRingDist);
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let matrix2 = sampler.sample_uniform(&params, 20, 5, DistType::FinRingDist);

        let sampler2 = DCRTPolyUniformSampler::new();
        let matrix3 = sampler2.sample_uniform(&params, 5, 12, DistType::FinRingDist);
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix arithmetic
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    #[test]
    fn test_gaussian_dist() {
        let params = DCRTPolyParams::default();

        // Test GaussianDist
        let sampler = DCRTPolyUniformSampler::new();
        let matrix1 =
            sampler.sample_uniform(&params, 20, 5, DistType::GaussDist { sigma: 4.57825 });
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let matrix2 =
            sampler.sample_uniform(&params, 20, 5, DistType::GaussDist { sigma: 4.57825 });

        let sampler2 = DCRTPolyUniformSampler::new();
        let matrix3 = sampler2.sample_uniform(&params, 5, 12, DistType::FinRingDist);
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix arithmetic
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    #[test]
    fn test_bit_dist() {
        let params = DCRTPolyParams::default();

        // Test BitDist
        let sampler = DCRTPolyUniformSampler::new();
        let matrix1 = sampler.sample_uniform(&params, 20, 5, DistType::BitDist);
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);
        // [TODO] Test the norm of each coefficient of polynomials in the matrix.

        let matrix2 = sampler.sample_uniform(&params, 20, 5, DistType::BitDist);

        let sampler2 = DCRTPolyUniformSampler::new();
        let matrix3 = sampler2.sample_uniform(&params, 5, 12, DistType::FinRingDist);
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix arithmetic
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }
}
