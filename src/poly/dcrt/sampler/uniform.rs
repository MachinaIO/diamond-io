#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
        sampler::{DistType, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use openfhe::ffi;

pub struct DCRTPolyUniformSampler {}

impl Default for DCRTPolyUniformSampler {
    fn default() -> Self {
        Self::new()
    }
}

impl DCRTPolyUniformSampler {
    pub fn new() -> Self {
        Self {}
    }

    pub fn sample_poly(&self, params: &DCRTPolyParams, dist: &DistType) -> DCRTPoly {
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
        DCRTPoly::new(sampled_poly)
    }
}

impl PolyUniformSampler for DCRTPolyUniformSampler {
    type M = DCRTPolyMatrix;

    fn sample_uniform(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        rows: usize,
        columns: usize,
        dist: DistType,
    ) -> Self::M {
        let c: Vec<Vec<DCRTPoly>> = parallel_iter!(0..rows)
            .map(|_| {
                parallel_iter!(0..columns)
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

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    #[cfg(not(feature = "parallel"))]
    use itertools::Itertools;
    #[cfg(not(feature = "parallel"))]
    use num_bigint::BigUint;
    #[cfg(not(feature = "parallel"))]
    use proptest::prelude::*;
    #[cfg(not(feature = "parallel"))]
    use std::sync::Arc;

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

    #[cfg(not(feature = "parallel"))]
    proptest::proptest! {
        #![proptest_config(ProptestConfig::with_cases(10))]

        #[test]
        fn test_bitdecomposition_uniform_sampler_ring(
            rows in 1usize..10usize,
            columns in 1usize..10usize,
        ) {
            let params = DCRTPolyParams::default();
            println!("Testing with: rows={}, columns={}, dist={:?}", rows, columns,  DistType::FinRingDist);
            let s = DCRTPolyUniformSampler::new();
            let matrix = s.sample_uniform(&params,rows, columns, DistType::FinRingDist);
            assert!(rows > 0 && columns > 0, "Invalid dimensions");
            let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, rows);
            let decomposed = matrix.decompose();
            let expected_matrix = gadget_matrix * decomposed;
            assert_eq!(matrix, expected_matrix);
        }

        #[test]
        fn test_bitdecomposition_uniform_sampler_bit(
            rows in 1usize..10usize,
            columns in 1usize..10usize,
        ) {
            let params = DCRTPolyParams::default();
            println!("Testing with: rows={}, columns={}, dist={:?}", rows, columns,  DistType::BitDist);
            let s = DCRTPolyUniformSampler::new();
            let matrix = s.sample_uniform(&params,rows, columns, DistType::BitDist);
            assert!(rows > 0 && columns > 0, "Invalid dimensions");
            let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, rows);
            let decomposed = matrix.decompose();
            let expected_matrix = gadget_matrix * decomposed;
            assert_eq!(matrix, expected_matrix);
        }

        #[test]
        fn test_bitdecomposition_uniform_sampler_gaussian(
            rows in 1usize..10usize,
            columns in 1usize..10usize,
            sigma in 1.0f64..10.0f64,
        ) {
            let params = DCRTPolyParams::default();
            println!("Testing with: rows={}, columns={}, dist={:?}", rows, columns, DistType::GaussDist { sigma });
            let s = DCRTPolyUniformSampler::new();
            let matrix = s.sample_uniform(&params,rows, columns, DistType::GaussDist { sigma });
            assert!(rows > 0 && columns > 0, "Invalid dimensions");
            let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, rows);
            let decomposed = matrix.decompose();
            let expected_matrix = gadget_matrix * decomposed;
            assert_eq!(matrix, expected_matrix);
        }

        #[test]
        fn test_modulus_switch_uniform_sampler_ring(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
        ) {
            let params = DCRTPolyParams::default();
            println!("Testing with: rows={}, columns={}, dist={:?}", rows, columns,  DistType::FinRingDist);
            let s = DCRTPolyUniformSampler::new();
            let matrix = s.sample_uniform(&params,rows, columns, DistType::FinRingDist);
            assert!(rows > 0 && columns > 0, "Invalid dimensions");

            let new_modulus = Arc::new(BigUint::from(2u32));
            let switched = matrix.modulus_switch(&new_modulus);
            let mut expected_matrix_vec = vec![];
            for i in 0..rows {
                let mut row = vec![];
                for j in 0..columns {
                    let poly = matrix.entry(i, j);
                    let coeffs = poly.coeffs().iter().map(|coeff| coeff.modulus_switch(new_modulus.clone())).collect_vec();
                    let new_poly = DCRTPoly::from_coeffs(&params, &coeffs);
                    row.push(new_poly);
                }
                expected_matrix_vec.push(row);
            }
            let expected = DCRTPolyMatrix::from_poly_vec(&params, expected_matrix_vec);
            assert_eq!(switched, expected);
        }

        #[test]
        fn test_modulus_switch_uniform_sampler_bit(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
        ) {
            let params = DCRTPolyParams::default();
            println!("Testing with: rows={}, columns={}, dist={:?}", rows, columns,  DistType::BitDist);
            let s = DCRTPolyUniformSampler::new();
            let matrix = s.sample_uniform(&params,rows, columns, DistType::BitDist);
            assert!(rows > 0 && columns > 0, "Invalid dimensions");

            let new_modulus = Arc::new(BigUint::from(2u32));
            let switched = matrix.modulus_switch(&new_modulus);
            let mut expected_matrix_vec = vec![];
            for i in 0..rows {
                let mut row = vec![];
                for j in 0..columns {
                    let poly = matrix.entry(i, j);
                    let coeffs = poly.coeffs().iter().map(|coeff| coeff.modulus_switch(new_modulus.clone())).collect_vec();
                    let new_poly = DCRTPoly::from_coeffs(&params, &coeffs);
                    row.push(new_poly);
                }
                expected_matrix_vec.push(row);
            }
            let expected = DCRTPolyMatrix::from_poly_vec(&params, expected_matrix_vec);
            assert_eq!(switched, expected);
        }

        #[test]
        fn test_modulus_switch_uniform_sampler_gaussian(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
            sigma in 1.0f64..10.0f64,
        ) {
            let params = DCRTPolyParams::default();
            println!("Testing with: rows={}, columns={}, dist={:?}", rows, columns,  DistType::GaussDist { sigma });
            let s = DCRTPolyUniformSampler::new();
            let matrix = s.sample_uniform(&params,rows, columns, DistType::GaussDist { sigma });
            assert!(rows > 0 && columns > 0, "Invalid dimensions");

            let new_params = DCRTPolyParams::new(4, 15, 51);
            let new_modulus = Arc::new(BigUint::from(2u32));
            let switched = matrix.modulus_switch(&new_modulus);
            let mut expected_matrix_vec = vec![];
            for i in 0..rows {
                let mut row = vec![];
                for j in 0..columns {
                    let poly = matrix.entry(i, j);
                    let coeffs = poly.coeffs().iter().map(|coeff| coeff.modulus_switch(new_modulus.clone())).collect_vec();
                    let new_poly = DCRTPoly::from_coeffs(&params, &coeffs);
                    row.push(new_poly);
                }
                expected_matrix_vec.push(row);
            }
            let expected = DCRTPolyMatrix::from_poly_vec(&new_params, expected_matrix_vec);
            assert_eq!(switched, expected);
        }
    }
}
