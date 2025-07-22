use super::{element::FinRingElem, params::DCRTPolyParams};
use crate::{
    impl_binop_with_refs, parallel_iter,
    poly::{element::PolyElem, Poly, PolyParams},
    utils::chunk_size_for,
};
use num_bigint::BigUint;
use num_traits::Zero;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPoly as DCRTPolyCxx},
    parse_coefficients_bytes,
};
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};

#[derive(Clone, Debug)]
pub struct DCRTPoly {
    ptr_poly: Arc<UniquePtr<DCRTPolyCxx>>,
}

// SAFETY: DCRTPoly is plain old data and is shared across threads in C++ OpenFHE as well.
unsafe impl Send for DCRTPoly {}
unsafe impl Sync for DCRTPoly {}

impl DCRTPoly {
    pub fn new(ptr_poly: UniquePtr<DCRTPolyCxx>) -> Self {
        Self { ptr_poly: ptr_poly.into() }
    }

    pub fn get_poly(&self) -> &UniquePtr<DCRTPolyCxx> {
        &self.ptr_poly
    }

    pub fn modulus_switch(
        &self,
        params: &DCRTPolyParams,
        new_modulus: <DCRTPolyParams as PolyParams>::Modulus,
    ) -> Self {
        debug_assert!(new_modulus < params.modulus());
        let coeffs = self.coeffs();
        let new_coeffs = coeffs
            .par_iter()
            .map(|coeff| coeff.modulus_switch(new_modulus.clone()))
            .collect::<Vec<FinRingElem>>();
        DCRTPoly::from_coeffs(params, &new_coeffs)
    }

    fn poly_gen_from_vec(params: &DCRTPolyParams, values: Vec<String>) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromVec(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            &values,
        ))
    }

    fn poly_gen_from_const(params: &DCRTPolyParams, value: String) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.ring_dimension(),
            params.crt_depth(),
            params.crt_bits(),
            &value,
        ))
    }
}

impl Poly for DCRTPoly {
    type Elem = FinRingElem;
    type Params = DCRTPolyParams;

    fn coeffs(&self) -> Vec<Self::Elem> {
        let poly_encoding = self.ptr_poly.GetCoefficientsBytes();
        let parsed_values = parse_coefficients_bytes(&poly_encoding);
        let coeffs = parsed_values.coefficients;
        let modulus = parsed_values.modulus;
        parallel_iter!(coeffs).map(|s| FinRingElem::new(s, Arc::new(modulus.clone()))).collect()
    }

    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Self {
        let mut coeffs_cxx = Vec::with_capacity(coeffs.len());
        for coeff in coeffs {
            debug_assert_eq!(coeff.modulus(), params.modulus().as_ref());
            coeffs_cxx.push(coeff.value().to_string());
        }
        Self::poly_gen_from_vec(params, coeffs_cxx)
    }

    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self {
        Self::poly_gen_from_const(params, constant.value().to_string())
    }

    fn from_decomposed(params: &DCRTPolyParams, decomposed: &[Self]) -> Self {
        let mut reconstructed = Self::const_zero(params);
        for (i, bit_poly) in decomposed.iter().enumerate() {
            let power_of_two = BigUint::from(2u32).pow(i as u32);
            let const_poly_power_of_two =
                Self::from_const(params, &FinRingElem::new(power_of_two, params.modulus()));
            reconstructed += bit_poly * &const_poly_power_of_two;
        }
        reconstructed
    }

    fn const_zero(params: &Self::Params) -> Self {
        Self::poly_gen_from_const(params, BigUint::ZERO.to_string())
    }

    fn const_one(params: &Self::Params) -> Self {
        Self::poly_gen_from_const(params, BigUint::from(1u32).to_string())
    }

    fn const_minus_one(params: &Self::Params) -> Self {
        Self::poly_gen_from_const(
            params,
            (params.modulus().as_ref() - BigUint::from(1u32)).to_string(),
        )
    }

    fn const_power_of_base(params: &Self::Params, k: usize) -> Self {
        let base = 1u32 << params.base_bits();
        Self::poly_gen_from_const(params, BigUint::from(base).pow(k as u32).to_string())
    }

    fn const_max(params: &Self::Params) -> Self {
        let coeffs = vec![FinRingElem::max_q(&params.modulus()); params.ring_dimension() as usize];
        Self::from_coeffs(params, &coeffs)
    }

    fn const_int(params: &Self::Params, int: usize) -> Self {
        Self::poly_gen_from_const(params, BigUint::from(int).to_string())
    }

    /// Encode `int` in little-endian bit order
    fn from_const_int_lsb(params: &Self::Params, int: usize) -> Self {
        let n = params.ring_dimension() as usize;
        let q = params.modulus();
        let one = FinRingElem::one(&q);
        let zero = FinRingElem::zero(&q);

        let coeffs: Vec<FinRingElem> =
            (0..n).map(|i| if (int >> i) & 1 == 1 { one.clone() } else { zero.clone() }).collect();

        Self::from_coeffs(params, &coeffs)
    }

    /// Decompose a polynomial of form b_0 + b_1 * x + b_2 * x^2 + ... + b_{n-1} * x^{n-1}
    /// where b_{j, h} is the h-th digit of the j-th coefficient of the polynomial.
    /// Return a vector of polynomials, where the h-th polynomial is defined as
    /// b_{0, h} + b_{1, h} * x + b_{2, h} * x^2 + ... + b_{n-1, h} * x^{n-1}.
    fn decompose_base(&self, params: &Self::Params) -> Vec<Self> {
        let coeffs = self.coeffs();
        let log_q = params.modulus_bits();
        let base_bits = params.base_bits() as usize;

        // Calculate the number of digits needed in the decomposition
        let num_digits = params.modulus_digits();

        // Create a mask for extracting the base_bits bits
        let base_mask = (BigUint::from(1u32) << base_bits) - BigUint::from(1u32);

        // Directly decompose into base_bits digits
        parallel_iter!(0..num_digits)
            .map(|digit_idx| {
                // Compute the shift amount for this digit
                let shift_amount = digit_idx * base_bits;

                // Extract the digit values for all coefficients
                let digit_values = coeffs
                    .par_iter()
                    .map(|coeff| {
                        if shift_amount >= log_q {
                            BigUint::from(0u32) // Handle the case where shift exceeds modulus bits
                        } else {
                            (coeff.value() >> shift_amount) & &base_mask
                        }
                    })
                    .collect::<Vec<_>>();

                // Create a polynomial from these digit values
                let poly_from_digits = DCRTPoly::from_coeffs(
                    params,
                    &digit_values
                        .par_iter()
                        .map(|value| FinRingElem::new(value.clone(), params.modulus()))
                        .collect::<Vec<_>>(),
                );

                poly_from_digits
            })
            .collect()
    }

    /// Create a polynomial from a compact byte representation based on `to_compact_bytes` encoding
    #[inline]
    fn from_compact_bytes(params: &Self::Params, bytes: &[u8]) -> Self {
        let ring_dimension = params.ring_dimension() as usize;
        let modulus = params.modulus();

        // Parse header.
        let max_byte_size = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;
        let bit_vector_byte_size = ring_dimension.div_ceil(8);
        let bit_vector = &bytes[4..4 + bit_vector_byte_size];
        let coeffs_base_offset = 4 + bit_vector_byte_size;
        // Remaining bytes contain coefficient values
        let coeffs: Vec<FinRingElem> = reconstruct_coeffs_chunked(
            bytes,
            ring_dimension,
            max_byte_size,
            bit_vector,
            coeffs_base_offset,
            &modulus,
            chunk_size_for(ring_dimension),
        );

        Self::from_coeffs(params, &coeffs)
    }

    /// Convert the polynomial to a compact byte representation
    /// The returned bytes vector is encoded as follows:
    /// 1. The first four bytes contain the `max_byte_size`, namely the maximum byte size of any
    ///    coefficient in the poly
    /// 2. The next `ceil(n/8)` bytes contain a bit vector, where each bit indicates if the
    ///    corresponding coefficient is negative (> `q_half`) and `n` is the ring dimension
    /// 3. The remaining `n * max_byte_size` contain the coefficient values
    #[inline]
    fn to_compact_bytes(&self) -> Vec<u8> {
        let coeffs = self.coeffs();
        let ring_dimension = coeffs.len();
        let processed_coeffs: Vec<(bool, Vec<u8>)> =
            process_coeffs_chunked(&coeffs, chunk_size_for(ring_dimension));

        build_compact_bytes(processed_coeffs, ring_dimension)
    }

    /// Recover bits from a polynomial using decision thresholds q/4 and 3q/4
    fn extract_bits_with_threshold(&self, params: &Self::Params) -> Vec<bool> {
        let modulus = params.modulus();
        let half_q = FinRingElem::half_q(&modulus); // q/2
        let quarter_q = half_q.value() >> 1; // q/4
        let three_quarter_q = &quarter_q * 3u32; // 3q/4

        self.coeffs()
            .iter()
            .map(|coeff| coeff.value())
            .map(|coeff| coeff >= &quarter_q && coeff < &three_quarter_q)
            .collect()
    }

    fn to_bool_vec(&self) -> Vec<bool> {
        self.coeffs()
            .into_iter()
            .map(|c| {
                let v = c.value();
                if v == &BigUint::from(0u32) {
                    false
                } else if v == &BigUint::from(1u32) {
                    true
                } else {
                    panic!("Coefficient is not 0 or 1: {v}");
                }
            })
            .collect()
    }

    fn to_const_int(&self) -> usize {
        let mut sum = 0;
        for (i, c) in self.coeffs_digits().into_iter().enumerate() {
            sum += 2_u32.pow(i as u32) * c;
        }
        sum as usize
    }

    fn from_bool_vec(params: &Self::Params, coeffs: &[bool]) -> Self {
        let coeffs: Vec<_> = coeffs
            .into_iter()
            .map(|i| FinRingElem::constant(&params.modulus(), *i as u64))
            .collect();
        DCRTPoly::from_coeffs(&params, &coeffs)
    }
}

impl PartialEq for DCRTPoly {
    fn eq(&self, other: &Self) -> bool {
        if self.ptr_poly.is_null() || other.ptr_poly.is_null() {
            return false;
        }
        self.ptr_poly.IsEqual(&other.ptr_poly)
    }
}

impl Eq for DCRTPoly {}

impl_binop_with_refs!(DCRTPoly => Add::add(self, rhs: &DCRTPoly) -> DCRTPoly {
    DCRTPoly::new(ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly))
});

impl_binop_with_refs!(DCRTPoly => Mul::mul(self, rhs: &DCRTPoly) -> DCRTPoly {
    DCRTPoly::new(ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly))
});

impl_binop_with_refs!(DCRTPoly => Sub::sub(self, rhs: &DCRTPoly) -> DCRTPoly {
    self + -rhs
});

impl Neg for DCRTPoly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        -&self
    }
}

impl Neg for &DCRTPoly {
    type Output = DCRTPoly;

    fn neg(self) -> Self::Output {
        DCRTPoly::new(self.ptr_poly.Negate())
    }
}

impl AddAssign for DCRTPoly {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl AddAssign<&DCRTPoly> for DCRTPoly {
    fn add_assign(&mut self, rhs: &Self) {
        // TODO: Expose `operator+=` in ffi.
        *self = &*self + rhs;
    }
}

impl MulAssign for DCRTPoly {
    fn mul_assign(&mut self, rhs: Self) {
        *self *= &rhs;
    }
}

impl MulAssign<&DCRTPoly> for DCRTPoly {
    fn mul_assign(&mut self, rhs: &Self) {
        // TODO: Expose `operator*=` in ffi.
        *self = &*self * rhs;
    }
}

impl SubAssign for DCRTPoly {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl SubAssign<&DCRTPoly> for DCRTPoly {
    fn sub_assign(&mut self, rhs: &Self) {
        *self += -rhs;
    }
}

// ==== Compact bytes for poly utils ====

fn build_compact_bytes(processed_coeffs: Vec<(bool, Vec<u8>)>, ring_dimension: usize) -> Vec<u8> {
    let bit_vector_byte_size = ring_dimension.div_ceil(8);
    let mut bit_vector = vec![0u8; bit_vector_byte_size];
    let mut max_byte_size = 1;

    // First pass: build bit vector and find max_byte_size.
    for (i, (is_negative, value_bytes)) in processed_coeffs.iter().enumerate() {
        if *is_negative {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            bit_vector[byte_idx] |= 1 << bit_idx;
        }
        max_byte_size = std::cmp::max(max_byte_size, value_bytes.len());
    }

    let total_byte_size = 4 + bit_vector_byte_size + (ring_dimension * max_byte_size);
    let mut result = vec![0u8; total_byte_size];

    // Write header.
    result[0..4].copy_from_slice(&(max_byte_size as u32).to_le_bytes());
    result[4..4 + bit_vector_byte_size].copy_from_slice(&bit_vector);

    // Write coefficient values.
    let coeffs_base_offset = 4 + bit_vector_byte_size;
    for (i, (_, value_bytes)) in processed_coeffs.iter().enumerate() {
        if !value_bytes.is_empty() {
            let start_pos = coeffs_base_offset + (i * max_byte_size);
            result[start_pos..start_pos + value_bytes.len()].copy_from_slice(value_bytes);
        }
    }

    result
}

fn reconstruct_coeffs_chunked(
    bytes: &[u8],
    ring_dimension: usize,
    max_byte_size: usize,
    bit_vector: &[u8],
    coeffs_base_offset: usize,
    modulus: &std::sync::Arc<BigUint>,
    chunk_size: usize,
) -> Vec<FinRingElem> {
    (0..ring_dimension)
        .into_par_iter()
        .chunks(chunk_size)
        .flat_map(|chunk| {
            chunk
                .into_iter()
                .map(|i| {
                    reconstruct_single_coeff(
                        i,
                        bytes,
                        max_byte_size,
                        bit_vector,
                        coeffs_base_offset,
                        modulus,
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

#[inline(always)]
fn reconstruct_single_coeff(
    i: usize,
    bytes: &[u8],
    max_byte_size: usize,
    bit_vector: &[u8],
    coeffs_base_offset: usize,
    modulus: &std::sync::Arc<BigUint>,
) -> FinRingElem {
    let start = coeffs_base_offset + (i * max_byte_size);

    // Find actual length by trimming trailing zeros.
    let mut value_len = max_byte_size;
    let value_bytes = &bytes[start..start + max_byte_size];
    while value_len > 0 && value_bytes[value_len - 1] == 0 {
        value_len -= 1;
    }

    // Only parse non-zero bytes.
    let value = if value_len == 0 {
        BigUint::ZERO
    } else {
        BigUint::from_bytes_le(&value_bytes[..value_len])
    };

    let byte_idx = i / 8;
    let bit_idx = i % 8;
    let is_negative = (bit_vector[byte_idx] & (1 << bit_idx)) != 0;

    // Convert back from centered representation.
    let final_value =
        if is_negative && !value.is_zero() { modulus.as_ref() - &value } else { value };

    FinRingElem::new(final_value, modulus.clone())
}

fn process_coeffs_chunked(coeffs: &[FinRingElem], chunk_size: usize) -> Vec<(bool, Vec<u8>)> {
    coeffs
        .par_chunks(chunk_size)
        .flat_map(|chunk| chunk.iter().map(|coeff| process_single_coeff(coeff)).collect::<Vec<_>>())
        .collect()
}

#[inline(always)]
fn process_single_coeff(coeff: &FinRingElem) -> (bool, Vec<u8>) {
    let coeff_val = coeff.value();
    let modulus_big = coeff.modulus();
    let q_half = modulus_big >> 1;

    if coeff_val > &q_half {
        let centered_value = coeff.modulus() - coeff_val;
        let value_bytes =
            if centered_value.is_zero() { Vec::new() } else { centered_value.to_bytes_le() };
        (true, value_bytes)
    } else if coeff_val.is_zero() {
        (false, Vec::new())
    } else {
        (false, coeff_val.to_bytes_le())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::DCRTPolyUniformSampler,
        sampler::{DistType, PolyUniformSampler},
        PolyParams,
    };
    use rand::prelude::*;
    use std::path::Path;

    #[test]
    fn test_const_int_roundtrip() {
        let mut rng = rand::rng();
        let params = DCRTPolyParams::default();

        for _ in 0..10 {
            let value = rng.random_range(0..(2_i32.pow(params.ring_dimension() - 1) as usize));
            let lsb_poly = DCRTPoly::from_const_int_lsb(&params, value);
            let poly = DCRTPoly::const_int(&params, value);
            let back = poly.to_const_int();
            let back_from_lsb = lsb_poly.to_const_int();
            assert_eq!(value, back);
            assert_eq!(value, back_from_lsb);
        }
    }

    #[test]
    fn test_dcrtpoly_coeffs() {
        let mut rng = rand::rng();
        /*
        todo: if x=0, n=1: libc++abi: terminating due to uncaught exception of type lbcrypto::OpenFHEException: /Users/piapark/Documents/GitHub/openfhe-development/src/core/include/math/nbtheory.h:l.156:ReverseBits(): msbb value not handled:0
        todo: if x=1, n=2: value mismatch from_coeffs & coeffs
        */
        let x = rng.random_range(12..20);
        let size = rng.random_range(1..20);
        let n = 2_i32.pow(x) as u32;
        let params = DCRTPolyParams::new(n, size, 51, 2);
        let q = params.modulus();
        let mut coeffs: Vec<FinRingElem> = Vec::new();
        for _ in 0..n {
            let value = rng.random_range(0..10000);
            coeffs.push(FinRingElem::new(value, q.clone()));
        }
        let poly = DCRTPoly::from_coeffs(&params, &coeffs);
        let extracted_coeffs = poly.coeffs();
        assert_eq!(coeffs, extracted_coeffs);
    }

    #[test]
    fn test_dcrtpoly_arithmetic() {
        let params = DCRTPolyParams::default();
        let q = params.modulus();

        // todo: replace value and modulus from param
        let coeffs1 = [
            FinRingElem::new(100u32, q.clone()),
            FinRingElem::new(200u32, q.clone()),
            FinRingElem::new(300u32, q.clone()),
            FinRingElem::new(400u32, q.clone()),
        ];
        let coeffs2 = [
            FinRingElem::new(500u32, q.clone()),
            FinRingElem::new(600u32, q.clone()),
            FinRingElem::new(700u32, q.clone()),
            FinRingElem::new(800u32, q.clone()),
        ];

        // 3. Create polynomials from those coefficients.
        let poly1 = DCRTPoly::from_coeffs(&params, &coeffs1);
        let poly2 = DCRTPoly::from_coeffs(&params, &coeffs2);

        // 4. Test addition.
        let sum = poly1.clone() + poly2.clone();

        // 5. Test multiplication.
        let product = &poly1 * &poly2;

        // 6. Test negation / subtraction.
        let neg_poly2 = poly2.clone().neg();
        let difference = poly1.clone() - poly2.clone();

        let mut poly_add_assign = poly1.clone();
        poly_add_assign += poly2.clone();

        let mut poly_mul_assign = poly1.clone();
        poly_mul_assign *= poly2.clone();

        // 8. Make some assertions
        assert!(sum != poly1, "Sum should differ from original poly1");
        assert!(neg_poly2 != poly2, "Negated polynomial should differ from original");
        assert_eq!(difference + poly2, poly1, "p1 - p2 + p2 should be p1");

        assert_eq!(poly_add_assign, sum, "+= result should match separate +");
        assert_eq!(poly_mul_assign, product, "*= result should match separate *");

        // 9. Test from_const / const_zero / const_one
        let const_poly = DCRTPoly::from_const(&params, &FinRingElem::new(123, q.clone()));
        assert_eq!(
            const_poly,
            DCRTPoly::from_coeffs(&params, &[FinRingElem::new(123, q.clone()); 1]),
            "from_const should produce a polynomial with all coeffs = 123"
        );
        let zero_poly = DCRTPoly::const_zero(&params);
        assert_eq!(
            zero_poly,
            DCRTPoly::from_coeffs(&params, &[FinRingElem::new(0, q.clone()); 1]),
            "const_zero should produce a polynomial with all coeffs = 0"
        );

        let one_poly = DCRTPoly::const_one(&params);
        assert_eq!(
            one_poly,
            DCRTPoly::from_coeffs(&params, &[FinRingElem::new(1, q); 1]),
            "one_poly should produce a polynomial with all coeffs = 1"
        );
    }

    #[test]
    fn test_dcrtpoly_decompose() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();
        let poly = sampler.sample_poly(&params, &DistType::FinRingDist);
        let decomposed = poly.decompose_base(&params);
        assert_eq!(decomposed.len(), { params.modulus_digits() });
    }

    #[test]
    fn test_dcrtpoly_to_compact_bytes_bit_dist() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();
        let poly = sampler.sample_poly(&params, &DistType::BitDist);
        let bytes = poly.to_compact_bytes();

        let ring_dimension = params.ring_dimension() as usize;

        // First four bytes contain `max_byte_size` (maximum byte size of any coefficient)
        // Since we're using BitDist, we expect max_byte_size to be equal to 1
        let max_byte_size = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;
        assert_eq!(max_byte_size, 1, "Max byte size size should be 1 for BitDist");

        // Next ceil(n/8) bytes are the bit vector (1 bit per coefficient)
        let bit_vector_byte_size = ring_dimension.div_ceil(8);

        // Expected total size
        // 4 bytes for `max_byte_size` + bit_vector_byte_size + (ring_dimension * max_byte_size)
        let expected_total_size = 4 + bit_vector_byte_size + (ring_dimension * max_byte_size);
        assert_eq!(bytes.len(), expected_total_size, "Incorrect total byte size");

        // Check that the structure is as expected
        // Verify bit vector section exists
        let bit_vector = &bytes[4..4 + bit_vector_byte_size];
        assert_eq!(bit_vector.len(), bit_vector_byte_size, "Bit vector size is incorrect");

        // Verify coefficient values section exists
        let coeffs_section = &bytes[4 + bit_vector_byte_size..];
        assert_eq!(
            coeffs_section.len(),
            ring_dimension * max_byte_size,
            "Coefficient section size is incorrect"
        );

        // Since we're using BitDist, each coefficient should be either 0 or 1
        // This means each byte in the coefficient section should be 0 or 1
        for (i, &coeff_byte) in coeffs_section.iter().enumerate() {
            assert!(
                coeff_byte == 0 || coeff_byte == 1,
                "Coefficient at position {} should be 0 or 1, got {}",
                i,
                coeff_byte
            );
        }
    }

    #[test]
    fn test_dcrtpoly_to_compact_bytes_uniform_dist() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();
        let poly = sampler.sample_poly(&params, &DistType::FinRingDist);
        let bytes = poly.to_compact_bytes();
        let modulus = params.modulus();
        let modulus_byte_size = modulus.to_bytes_le().len();

        let ring_dimension = params.ring_dimension() as usize;

        // First four bytes contain `max_byte_size` (maximum byte size of any coefficient)
        // Since we're using BitDist, we expect max_byte_size to be equal less or equal the byte
        // size of the modulus
        let max_byte_size = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;
        assert!(
            max_byte_size <= modulus_byte_size,
            "Max byte size should be less than or equal to modulus byte size"
        );

        // Next ceil(n/8) bytes are the bit vector (1 bit per coefficient)
        let bit_vector_byte_size = ring_dimension.div_ceil(8);

        // Expected total size
        // 4 bytes for `max_byte_size` + bit_vector_byte_size + (ring_dimension * max_byte_size)
        let expected_total_size = 4 + bit_vector_byte_size + (ring_dimension * max_byte_size);
        assert_eq!(bytes.len(), expected_total_size, "Incorrect total byte size");

        // Check that the structure is as expected
        // Verify bit vector section exists
        let bit_vector = &bytes[4..4 + bit_vector_byte_size];
        assert_eq!(bit_vector.len(), bit_vector_byte_size, "Bit vector size is incorrect");

        // Verify coefficient values section exists
        let coeffs_section = &bytes[4 + bit_vector_byte_size..];
        assert_eq!(
            coeffs_section.len(),
            ring_dimension * max_byte_size,
            "Coefficient section size is
        incorrect"
        );
    }

    #[test]
    fn test_dcrtpoly_from_compact_bytes() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Test with BitDist (binary coefficients)
        let original_poly = sampler.sample_poly(&params, &DistType::BitDist);
        let bytes = original_poly.to_compact_bytes();
        let reconstructed_poly = DCRTPoly::from_compact_bytes(&params, &bytes);

        // The original and reconstructed polynomials should be equal
        assert_eq!(
            original_poly, reconstructed_poly,
            "Reconstructed polynomial does not match original (BitDist)"
        );

        // Test with FinRingDist (random coefficients in the ring)
        let original_poly = sampler.sample_poly(&params, &DistType::FinRingDist);
        let bytes = original_poly.to_compact_bytes();
        let reconstructed_poly = DCRTPoly::from_compact_bytes(&params, &bytes);

        // The original and reconstructed polynomials should be equal
        assert_eq!(
            original_poly, reconstructed_poly,
            "Reconstructed polynomial does not match original (FinRingDist)"
        );

        // Test with GaussDist (Gaussian distribution)
        let original_poly = sampler.sample_poly(&params, &DistType::GaussDist { sigma: 3.2 });
        let bytes = original_poly.to_compact_bytes();
        let reconstructed_poly = DCRTPoly::from_compact_bytes(&params, &bytes);

        // The original and reconstructed polynomials should be equal
        assert_eq!(
            original_poly, reconstructed_poly,
            "Reconstructed polynomial does not match original (GaussDist)"
        );
    }

    #[tokio::test]
    async fn test_dcrtpoly_write_read_file() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();

        // Test with different distribution types
        let dists = [DistType::BitDist, DistType::FinRingDist, DistType::GaussDist { sigma: 3.0 }];

        // Create a temporary directory for testing
        let test_dir = Path::new("test_poly_write_read");
        if !test_dir.exists() {
            std::fs::create_dir(test_dir).unwrap();
        } else {
            // Clean it first to ensure no old files interfere
            std::fs::remove_dir_all(test_dir).unwrap();
            std::fs::create_dir(test_dir).unwrap();
        }

        for dist in dists {
            // Create a random polynomial
            let poly = sampler.sample_poly(&params, &dist);
            let poly_id = format!("test_poly_{:?}", dist);

            // Write the polynomial to a file
            poly.write_to_file(test_dir, &poly_id).await;

            // Read the polynomial back
            let read_poly = DCRTPoly::read_from_file(&params, test_dir, &poly_id);

            // Verify the polynomials are equal
            assert_eq!(poly, read_poly, "Read polynomial does not match original for {:?}", dist);
        }

        // Clean up
        std::fs::remove_dir_all(test_dir).unwrap();
    }
}
