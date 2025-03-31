use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, FinRingElem},
        sampler::{DistType, PolyHashSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use bitvec::prelude::*;
use digest::OutputSizeUser;
use num_bigint::BigUint;
use num_traits::Zero;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::{marker::PhantomData, ops::Range};

pub struct DCRTPolyHashSampler<H: OutputSizeUser + digest::Digest> {
    key: [u8; 32],
    _h: PhantomData<H>,
}

impl<H> DCRTPolyHashSampler<H>
where
    H: OutputSizeUser + digest::Digest,
{
    pub fn new(key: [u8; 32]) -> Self {
        Self { key, _h: PhantomData }
    }

    fn _set_key(&mut self, key: [u8; 32]) {
        self.key = key
    }

    fn _expose_key(&self) -> &[u8] {
        &self.key
    }
}

impl<H> PolyHashSampler<[u8; 32]> for DCRTPolyHashSampler<H>
where
    H: OutputSizeUser + digest::Digest + Clone + Send + Sync,
{
    type M = DCRTPolyMatrix;

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        tag: B,
        nrow: usize,
        ncol: usize,
        dist: DistType,
    ) -> DCRTPolyMatrix {
        let hash_output_size = <H as digest::Digest>::output_size() * 8;
        let n = params.ring_dimension() as usize;
        let q = params.modulus();
        let log_q = params.modulus_bits();
        let num_hash_fin_per_poly = (log_q * n).div_ceil(hash_output_size);
        let num_hash_bit_per_poly = n.div_ceil(hash_output_size);
        let mut new_matrix = DCRTPolyMatrix::new_empty(params, nrow, ncol);
        let mut hasher: H = H::new();
        hasher.update(self.key);
        hasher.update(tag.as_ref());
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
            match dist {
                DistType::FinRingDist => parallel_iter!(row_offsets)
                    .map(|i| {
                        parallel_iter!(col_offsets.clone())
                            .map(|j| {
                                let mut hasher = hasher.clone();
                                hasher.update(i.to_le_bytes());
                                hasher.update(j.to_le_bytes());
                                let mut local_bits = bitvec![u8, Lsb0;];
                                for hash_idx in 0..num_hash_fin_per_poly {
                                    let mut hasher = hasher.clone();
                                    hasher.update((hash_idx as u64).to_le_bytes());
                                    for &byte in hasher.finalize().iter() {
                                        for bit_index in 0..8 {
                                            local_bits.push((byte >> bit_index) & 1 != 0);
                                        }
                                    }
                                }
                                let local_bits = local_bits.split_at(log_q * n).0;
                                let coeffs = parallel_iter!(0..n)
                                    .map(|coeff_idx| {
                                        let bits =
                                            &local_bits[coeff_idx * log_q..(coeff_idx + 1) * log_q];
                                        let mut value = BigUint::zero();
                                        for bit in bits.iter() {
                                            value <<= 1;
                                            if *bit {
                                                value |= BigUint::from(1u32);
                                            }
                                        }
                                        FinRingElem::new(value, q.clone())
                                    })
                                    .collect::<Vec<_>>();
                                DCRTPoly::from_coeffs(params, &coeffs)
                            })
                            .collect()
                    })
                    .collect::<Vec<Vec<DCRTPoly>>>(),
                DistType::BitDist => parallel_iter!(row_offsets)
                    .map(|i| {
                        parallel_iter!(col_offsets.clone())
                            .map(|j| {
                                let mut hasher = hasher.clone();
                                hasher.update(i.to_le_bytes());
                                hasher.update(j.to_le_bytes());
                                let mut local_bits = bitvec![u8, Lsb0;];
                                for hash_idx in 0..num_hash_bit_per_poly {
                                    let mut hasher = hasher.clone();
                                    hasher.update((hash_idx as u64).to_le_bytes());
                                    for &byte in hasher.finalize().iter() {
                                        for bit_index in 0..8 {
                                            local_bits.push((byte >> bit_index) & 1 != 0);
                                        }
                                    }
                                }
                                let local_bits = local_bits.split_at(n).0;
                                let coeffs = parallel_iter!(0..n)
                                    .map(|coeff_idx| {
                                        FinRingElem::new(local_bits[coeff_idx] as u64, q.clone())
                                    })
                                    .collect::<Vec<_>>();
                                DCRTPoly::from_coeffs(params, &coeffs)
                            })
                            .collect::<Vec<DCRTPoly>>()
                    })
                    .collect::<Vec<Vec<DCRTPoly>>>(),
                _ => {
                    panic!("Unsupported distribution type")
                }
            }
        };
        new_matrix.replace_entries(0..nrow, 0..ncol, f);
        new_matrix
    }

    fn set_key(&mut self, key: [u8; 32]) {
        self._set_key(key)
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;
    #[cfg(not(feature = "parallel"))]
    use itertools::Itertools;
    use keccak_asm::Keccak256;
    #[cfg(not(feature = "parallel"))]
    use proptest::prelude::*;
    #[cfg(not(feature = "parallel"))]
    use std::sync::Arc;

    #[test]
    fn test_poly_hash_sampler() {
        let key = [0u8; 32];
        let params = DCRTPolyParams::default();
        let mut sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result = sampler.sample_hash(&params, tag, nrow, ncol, DistType::BitDist);
        // [TODO] Test the norm of each coefficient of polynomials in the matrix.

        let new_key = [1u8; 32];
        sampler.set_key(new_key);

        let matrix = matrix_result;
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }

    #[test]
    fn test_poly_hash_sampler_fin_ring_dist() {
        let key = [0u8; 32];
        let params = DCRTPolyParams::default();
        let mut sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result = sampler.sample_hash(&params, tag, nrow, ncol, DistType::FinRingDist);

        let new_key = [1u8; 32];
        sampler.set_key(new_key);

        let matrix = matrix_result;
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }

    #[cfg(not(feature = "parallel"))]
    proptest! {
        #![proptest_config(ProptestConfig::with_cases(10))]

        #[test]
        fn test_bitdecomposition_hash_sampler_ring(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
            key in any::<[u8; 32]>(),
            tag in any::<u64>(),
        ) {
            let params = DCRTPolyParams::default();
            let tag_bytes = tag.to_le_bytes();
            let sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
            let matrix = sampler.sample_hash(&params,tag_bytes, rows, columns, DistType::FinRingDist);
            let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, rows);
            let decomposed = matrix.decompose();
            let expected_matrix = gadget_matrix * decomposed;
            assert_eq!(matrix, expected_matrix);
        }

        #[test]
        fn test_bitdecomposition_hash_sampler_bit(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
            key in any::<[u8; 32]>(),
            tag in any::<u64>(),
        ) {
            let params = DCRTPolyParams::default();
            let tag_bytes = tag.to_le_bytes();
            let sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
            let matrix = sampler.sample_hash(&params,tag_bytes, rows, columns, DistType::BitDist);
            let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, rows);
            let decomposed = matrix.decompose();
            let expected_matrix = gadget_matrix * decomposed;
            assert_eq!(matrix, expected_matrix);
        }

        #[test]
        fn test_modulus_switch_hash_sampler_ring(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
            key in any::<[u8; 32]>(),
            tag in any::<u64>(),
        ) {
            let params = DCRTPolyParams::default();

            let tag_bytes = tag.to_le_bytes();
            let sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
            let matrix = sampler.sample_hash(&params,tag_bytes, rows, columns, DistType::FinRingDist);

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
        fn test_modulus_switch_hash_sampler_bit(
            rows in 1usize..5usize,
            columns in 1usize..5usize,
            key in any::<[u8; 32]>(),
            tag in any::<u64>(),
        ) {
            let params = DCRTPolyParams::default();

            let tag_bytes = tag.to_le_bytes();
            let sampler = DCRTPolyHashSampler::<Keccak256>::new(key);
            let matrix = sampler.sample_hash(&params,tag_bytes, rows, columns, DistType::BitDist);

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
    }
}
