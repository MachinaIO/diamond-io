use std::marker::PhantomData;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, FinRingElem},
    sampler::{DistType, PolyHashSampler},
    Poly, PolyMatrix, PolyParams,
};
use digest::OutputSizeUser;
use num_bigint::BigUint;

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

        // From field elements to nrow * ncol polynomials
        let total_poly = nrow * ncol;
        let mut offset = 0;
        let mut all_polys = Vec::with_capacity(total_poly);
        for _ in 0..total_poly {
            let coeffs = &ring_elems[offset..offset + n];
            offset += n;
            let poly = DCRTPoly::from_coeffs(params, coeffs);
            all_polys.push(poly);
        }

        // From polynomials to matrix such that the first row of the matrixcontains the first ncol
        // polynomials
        let mut matrix_inner = Vec::with_capacity(nrow);
        let mut poly_iter = all_polys.into_iter();
        for _ in 0..nrow {
            let row_polys: Vec<DCRTPoly> = poly_iter.by_ref().take(ncol).collect();
            matrix_inner.push(row_polys);
        }

        DCRTPolyMatrix::from_poly_vec(params, matrix_inner)
    }
}

impl<H> PolyHashSampler<[u8; 32]> for DCRTPolyHashSampler<H>
where
    H: OutputSizeUser + digest::Digest,
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
                // index = number of hashes to be performed = ceil(nrow * ncol * n * ceil(log2(q)) / hash_output_size)
                // field_elements = number of field elements sampled = (bits / ceil(log2(q)))
                let bit_length = params.modulus_bits();
                let index = (nrow * ncol * n * bit_length).div_ceil(hash_output_size);
                // bits = number of resulting bits from hashing ops = hash_output_size * index
                let mut bits = Vec::with_capacity(hash_output_size * index);
                let mut ring_elems = Vec::with_capacity((index * hash_output_size) / bit_length);
                for i in 0..(index) {
                    //  H ( key || tag || i )
                    let mut hasher = H::new();
                    // todo: we currently assuming index is less than u32
                    let min_i_type = std::mem::size_of_val(&i);
                    if min_i_type > std::mem::size_of::<u8>() {
                        let mut combined =
                            Vec::with_capacity(self.key.len() + tag.as_ref().len() + 32);
                        combined.extend_from_slice(&self.key);
                        combined.extend_from_slice(tag.as_ref());
                        combined.extend_from_slice(&i.to_be_bytes());
                        hasher.update(&combined);
                    } else {
                        let mut combined =
                            Vec::with_capacity(self.key.len() + tag.as_ref().len() + 1);
                        combined.extend_from_slice(&self.key);
                        combined.extend_from_slice(tag.as_ref());
                        combined.push(i as u8);
                        hasher.update(&combined);
                    }

                    for &byte in hasher.finalize().iter() {
                        for bit_index in 0..8 {
                            let bit = (byte >> bit_index) & 1;
                            bits.push(bit);
                        }
                    }
                }
                // From bits to field elements
                let mut offset = 0;
                for _ in 0..(bits.len() / bit_length) {
                    let value_bits = &bits[offset..offset + bit_length];
                    let value = BigUint::from_radix_be(value_bits, 2).unwrap();
                    offset += bit_length;
                    let fe = FinRingElem::new(value, q.clone());
                    ring_elems.push(fe);
                }
                ring_elems
            }
            DistType::BitDist => {
                // index = number of hashes to be performed = ceil(nrow * ncol * n * ceil(log2(q)) / hash_output_size)
                // field_elements = number of field elements sampled = (bits / ceil(log2(q)))
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
    use itertools::Itertools;
    use keccak_asm::Keccak256;
    use proptest::prelude::*;

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

            let new_params = DCRTPolyParams::new(4, 15, 51);
            let new_modulus = new_params.modulus();
            let switched = matrix.modulus_switch(&new_params);
            assert_eq!(switched.params().modulus(), new_params.modulus());
            let mut expected_matrix_vec = vec![];
            for i in 0..rows {
                let mut row = vec![];
                for j in 0..columns {
                    let poly = matrix.entry(i, j);
                    let coeffs = poly.coeffs().iter().map(|coeff| coeff.modulus_switch(new_modulus.clone())).collect_vec();
                    let new_poly = DCRTPoly::from_coeffs(&new_params, &coeffs);
                    row.push(new_poly);
                }
                expected_matrix_vec.push(row);
            }
            let expected = DCRTPolyMatrix::from_poly_vec(&new_params, expected_matrix_vec);
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

            let new_params = DCRTPolyParams::new(4, 15, 51);
            let new_modulus = new_params.modulus();
            let switched = matrix.modulus_switch(&new_params);
            assert_eq!(switched.params().modulus(), new_params.modulus());
            let mut expected_matrix_vec = vec![];
            for i in 0..rows {
                let mut row = vec![];
                for j in 0..columns {
                    let poly = matrix.entry(i, j);
                    let coeffs = poly.coeffs().iter().map(|coeff| coeff.modulus_switch(new_modulus.clone())).collect_vec();
                    let new_poly = DCRTPoly::from_coeffs(&new_params, &coeffs);
                    row.push(new_poly);
                }
                expected_matrix_vec.push(row);
            }
            let expected = DCRTPolyMatrix::from_poly_vec(&new_params, expected_matrix_vec);
            assert_eq!(switched, expected);
        }
    }
}
