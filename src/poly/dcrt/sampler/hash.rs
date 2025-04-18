use crate::{
    parallel_iter,
    poly::{
        Poly, PolyMatrix, PolyParams,
        dcrt::{DCRTPoly, DCRTPolyMatrix, FinRingElem},
        sampler::{DistType, PolyHashSampler},
    },
};
use bitvec::prelude::*;
use digest::OutputSizeUser;
use num_bigint::BigUint;
use num_traits::Zero;
use rayon::prelude::*;
use std::{marker::PhantomData, ops::Range};

pub struct DCRTPolyHashSampler<H: OutputSizeUser + digest::Digest> {
    _h: PhantomData<H>,
}

impl<H> PolyHashSampler<[u8; 32]> for DCRTPolyHashSampler<H>
where
    H: OutputSizeUser + digest::Digest + Clone + Send + Sync,
{
    type M = DCRTPolyMatrix;

    fn new() -> Self {
        Self { _h: PhantomData }
    }

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        params: &<<Self::M as PolyMatrix>::P as Poly>::Params,
        hash_key: [u8; 32],
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
        hasher.update(hash_key);
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;
    use keccak_asm::Keccak256;

    #[test]
    fn test_poly_hash_sampler() {
        let key = [0u8; 32];
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyHashSampler::<Keccak256>::new();
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result = sampler.sample_hash(&params, key, tag, nrow, ncol, DistType::BitDist);
        // [TODO] Test the norm of each coefficient of polynomials in the matrix.

        let matrix = matrix_result;
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }

    #[test]
    fn test_poly_hash_sampler_fin_ring_dist() {
        let key = [0u8; 32];
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyHashSampler::<Keccak256>::new();
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result =
            sampler.sample_hash(&params, key, tag, nrow, ncol, DistType::FinRingDist);

        let matrix = matrix_result;
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }
}
