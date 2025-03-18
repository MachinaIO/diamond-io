use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, FinRingElem},
        sampler::{DistType, PolyHashSampler},
        Poly, PolyMatrix, PolyParams,
    },
    utils::log_mem,
};
use bitvec::prelude::*;
use digest::OutputSizeUser;
use num_bigint::BigUint;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::marker::PhantomData;
use tracing::info;

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

    fn ring_elems_to_matrix(
        params: &DCRTPolyParams,
        ring_elems: Vec<FinRingElem>,
        nrow: usize,
        ncol: usize,
    ) -> DCRTPolyMatrix {
        let n = params.ring_dimension() as usize;
        // Check if we have enough field elements (not sure if this will ever happen)
        if ring_elems.len() < nrow * ncol * n {
            panic!("Not enough ring elements to sample hash")
        }

        info!("Begin ring_elems_to_matrix");
        log_mem();
        // From field elements to nrow * ncol polynomials
        DCRTPolyMatrix::from_poly_vec(
            params,
            parallel_iter!(0..nrow)
                .map(|row_idx| {
                    let row_offset = row_idx * ncol * n;
                    parallel_iter!(0..ncol)
                        .map(|col_idx| {
                            let offset = row_offset + col_idx * n;
                            let coeffs = &ring_elems[offset..offset + n];
                            DCRTPoly::from_coeffs(params, coeffs)
                        })
                        .collect()
                })
                .collect(),
        )
    }
}

impl<H> PolyHashSampler<[u8; 32]> for DCRTPolyHashSampler<H>
where
    H: OutputSizeUser + digest::Digest + Clone,
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

        let ring_elems = match dist {
            DistType::FinRingDist => {
                // index = number of hashes to be performed = ceil(nrow * ncol * n * ceil(log2(q)) /
                // hash_output_size) field_elements = number of field elements
                // sampled = (bits / ceil(log2(q)))
                let bit_length = params.modulus_bits();
                let index = (nrow * ncol * n * bit_length).div_ceil(hash_output_size);
                // bits = number of resulting bits from hashing ops = hash_output_size * index
                let mut bv = bitvec![u8, Lsb0;];
                let mut og_hasher: H = H::new();
                og_hasher.update(self.key);
                og_hasher.update(tag.as_ref());
                for i in 0..index {
                    let mut hasher = og_hasher.clone();
                    //  H ( key || tag || i )
                    hasher.update(i.to_le_bytes());
                    for &byte in hasher.finalize().iter() {
                        for bit_index in 0..8 {
                            bv.push((byte >> bit_index) & 1 != 0);
                        }
                    }
                }
                let num_chunks = bv.len() / bit_length;
                let ring_elems: Vec<FinRingElem> = parallel_iter!(0..num_chunks)
                    .map(|i| {
                        let start = i * bit_length;
                        let end = start + bit_length;
                        let value_bits = &bv[start..end];

                        let mut value = BigUint::default();
                        for bit in value_bits.iter() {
                            value <<= 1;
                            if *bit {
                                value |= BigUint::from(1u32);
                            }
                        }
                        FinRingElem::new(value, q.clone())
                    })
                    .collect();
                debug_assert_eq!(ring_elems.len(), (index * hash_output_size) / bit_length);
                ring_elems
            }
            DistType::BitDist => {
                // index = number of hashes to be performed = ceil(nrow * ncol * n * ceil(log2(q)) /
                // hash_output_size) field_elements = number of field elements
                // sampled = (bits / ceil(log2(q)))
                info!("Sampling bit distribution");
                log_mem();
                let index = (nrow * ncol * n).div_ceil(hash_output_size);
                let mut ring_elems = Vec::with_capacity(hash_output_size * index);
                for i in 0..index {
                    //  H ( key || tag || i )
                    let mut hasher = H::new();
                    let mut combined = Vec::with_capacity(self.key.len() + tag.as_ref().len() + 1);
                    combined.extend_from_slice(&self.key);
                    combined.extend_from_slice(tag.as_ref());
                    combined.push(i as u8);
                    hasher.update(&combined);
                    for &byte in hasher.finalize().iter() {
                        for bit_index in 0..8 {
                            info!("Sampling bit {}", bit_index);
                            log_mem();
                            let bit = (byte >> bit_index) & 1;
                            ring_elems.push(FinRingElem::new(bit as u64, q.clone()));
                        }
                    }
                }
                ring_elems
            }
            _ => {
                panic!("Unsupported distribution type")
            }
        };

        Self::ring_elems_to_matrix(params, ring_elems, nrow, ncol)
    }

    fn set_key(&mut self, key: [u8; 32]) {
        self._set_key(key)
    }
}

#[cfg(test)]
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
