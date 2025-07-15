use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, FinRingElem},
        sampler::{DistType, PolyHashSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use bitvec::prelude::*;
use digest::OutputSizeUser;
use num_bigint::BigUint;
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
        let out_sz = <H as digest::Digest>::output_size();
        let n = params.ring_dimension() as usize;
        let bytes_per_poly = (n + 7) / 8;
        let hashes_per_poly = (bytes_per_poly + out_sz - 1) / out_sz;
        let hash_output_size = <H as digest::Digest>::output_size() * 8;
        let q = params.modulus();
        let log_q = params.modulus_bits();
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
                                let mut buf = vec![0u8; hashes_per_poly * out_sz];
                                for blk in 0..hashes_per_poly {
                                    let mut h = hasher.clone(); // cheap 200‑B copy
                                    h.update((blk as u64).to_le_bytes()); // counter
                                    h.finalize_into(
                                        (&mut buf[blk * out_sz..(blk + 1) * out_sz]).into(),
                                    );
                                }

                                let bits = &buf[..bytes_per_poly];
                                let coeffs = if log_q <= 64 {
                                    (0..n)
                                        .map(|k| {
                                            let byte = bits[k >> 3]; // k / 8
                                            let bit = (byte >> (k & 7)) & 1;
                                            FinRingElem::new(bit as u64, q.clone())
                                        })
                                        .collect::<Vec<_>>()
                                } else {
                                    (0..n)
                                        .map(|k| {
                                            let byte = bits[k >> 3];
                                            let bit = (byte >> (k & 7)) & 1;
                                            FinRingElem::new(BigUint::from(bit), q.clone())
                                        })
                                        .collect::<Vec<_>>()
                                };
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
