use itertools::Itertools;
use rayon::prelude::*;

use super::{element::FinRingElem, params::DCRTPolyParams};
use crate::{
    impl_binop_with_refs, parallel_iter,
    poly::{element::PolyElem, Poly, PolyParams},
};
use num_bigint::BigUint;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPoly as DCRTPolyCxx},
    parse_coefficients_bytes,
};

use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
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

    /// Create a polynomial from a compact byte representation
    /// The bytes should have the format created by to_compact_bytes
    fn from_compact_bytes(params: &Self::Params, bytes: &[u8]) -> Self {
        let ring_dimension = params.ring_dimension() as usize;
        let modulus = params.modulus();

        // First four bytes contain the byte size per coefficient
        let max_byte_size = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;

        // Next n/8 bytes contain the bit vector
        let bit_vector_byte_size = ring_dimension.div_ceil(8);
        let bit_vector = &bytes[4..4 + bit_vector_byte_size];

        // Remaining bytes contain coefficient values
        let coeffs: Vec<FinRingElem> = parallel_iter!(0..ring_dimension)
            .map(|i| {
                let start = 4 + bit_vector_byte_size + (i * max_byte_size);
                let end = start + max_byte_size;
                let value_bytes = &bytes[start..end];

                let value = BigUint::from_bytes_le(value_bytes);

                let byte_idx = i / 8;
                let bit_idx = i % 8;
                let is_negative = (bit_vector[byte_idx] & (1 << bit_idx)) != 0;

                // Convert back from centered representation
                let final_value = if is_negative {
                    // If negative flag is set, compute q - value
                    modulus.as_ref() - &value
                } else {
                    // Otherwise, use value as is
                    value
                };

                FinRingElem::new(final_value, modulus.clone())
            })
            .collect();

        Self::from_coeffs(params, &coeffs)
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

    fn const_power_of_two(params: &Self::Params, k: usize) -> Self {
        Self::poly_gen_from_const(params, BigUint::from(2u32).pow(k as u32).to_string())
    }

    fn const_max(params: &Self::Params) -> Self {
        let coeffs = vec![FinRingElem::max_q(&params.modulus()); params.ring_dimension() as usize];
        Self::from_coeffs(params, &coeffs)
    }

    /// Decompose a polynomial of form b_0 + b_1 * x + b_2 * x^2 + ... + b_{n-1} * x^{n-1}
    /// where b_{j, h} is the h-th bit of the j-th coefficient of the polynomial.
    /// Return a vector of polynomials, where the h-th polynomial is defined as
    /// b_{0, h} + b_{1, h} * x + b_{2, h} * x^2 + ... + b_{n-1, h} * x^{n-1}.
    fn decompose_bits(&self, params: &Self::Params) -> Vec<Self> {
        let coeffs = self.coeffs();
        let bit_length = params.modulus_bits();
        parallel_iter!(0..bit_length)
            .map(|h| {
                DCRTPoly::from_coeffs(
                    params,
                    &coeffs
                        .iter()
                        .map(|j| {
                            FinRingElem::new(
                                (j.value() >> h) & BigUint::from(1u32),
                                params.modulus(),
                            )
                        })
                        .collect_vec(),
                )
            })
            .collect()
    }

    /// Convert the polynomial to a compact byte representation
    /// The returned bytes vector is encoded as follows:
    /// 1. The first four bytes contain the `max_byte_size`, namely the maximum byte size of any
    ///    coefficient in the poly
    /// 2. The next `ceil(n/8)` bytes contain a bit vector, where each bit indicates if the
    ///    corresponding coefficient is negative and `n` is the ring dimension
    /// 3. The remaining `n * max_byte_size` contain the coefficient values
    fn to_compact_bytes(&self) -> Vec<u8> {
        let modulus = self.ptr_poly.GetModulus();
        let modulus_big: BigUint = BigUint::from_str(&modulus).unwrap();
        let q_half = &modulus_big / 2u8;

        let coeffs = self.coeffs();
        let ring_dimension = coeffs.len();

        // Create a bit vector to store flags for negative coefficients
        let bit_vector_byte_size = ring_dimension.div_ceil(8);
        let mut bit_vector = vec![0u8; bit_vector_byte_size];

        let mut max_byte_size = 0;
        let mut processed_values = Vec::with_capacity(ring_dimension);

        // First pass: Process coefficients and calculate max_byte_size
        for (i, coeff) in coeffs.iter().enumerate() {
            // Center coefficients around 0
            let value = if coeff.value() > &q_half {
                let byte_idx = i / 8;
                let bit_idx = i % 8;
                bit_vector[byte_idx] |= 1 << bit_idx; // Set flag for negative coefficient
                &modulus_big - coeff.value() // Convert to absolute value: q - a
            } else if coeff.value() == &BigUint::ZERO {
                BigUint::ZERO
            } else {
                coeff.value().clone()
            };

            processed_values.push(value.clone());

            let value_bytes = value.to_bytes_le();
            max_byte_size = std::cmp::max(max_byte_size, value_bytes.len());
        }

        let total_byte_size = 4 + bit_vector_byte_size + (ring_dimension * max_byte_size);
        let mut result = vec![0u8; total_byte_size];

        // Store max_byte_size in the first four bytes (little-endian)
        let max_byte_size_bytes = (max_byte_size as u32).to_le_bytes();
        result[0..4].copy_from_slice(&max_byte_size_bytes);

        // Store bit vector
        result[4..4 + bit_vector_byte_size].copy_from_slice(&bit_vector);

        // Second pass: Store preprocessed coefficient values s.t. each coefficient is max_byte_size
        // bytes long
        for (i, value) in processed_values.iter().enumerate() {
            let value_bytes = value.to_bytes_le();
            let start_pos = 4 + bit_vector_byte_size + (i * max_byte_size);

            result[start_pos..start_pos + value_bytes.len()].copy_from_slice(&value_bytes);
        }

        result
    }

    /// Recover bits from a polynomial using decision thresholds q/4 and 3q/4
    fn extract_bits_with_threshold(&self, params: &Self::Params) -> Vec<bool> {
        let modulus = params.modulus();
        let quarter_q = modulus.as_ref() >> 2; // q/4
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
                    panic!("Coefficient is not 0 or 1: {}", v);
                }
            })
            .collect()
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

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::poly::{dcrt::DCRTPolyUniformSampler, sampler::DistType, PolyParams};
    use rand::prelude::*;

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
        let decomposed = poly.decompose_bits(&params);
        assert_eq!(decomposed.len(), params.modulus_bits());
    }

    #[test]
    fn test_dcrtpoly_to_compact_bytes() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();
        let poly = sampler.sample_poly(&params, &DistType::BitDist);
        let bytes = poly.to_compact_bytes();

        let ring_dimension = params.ring_dimension() as usize;

        // First four bytes contain `max_byte_size` (maximum byte size of any coefficient)
        // Since we're using BitDist, we expect max_byte_size to be equal to 1
        let max_byte_size = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
        assert_eq!(max_byte_size, 1, "Max byte size size should be 1 for BitDist");

        // Next ceil(n/8) bytes are the bit vector (1 bit per coefficient)
        let bit_vector_byte_size = ring_dimension.div_ceil(8);

        // Expected total size
        // 4 bytes for `max_byte_size` + bit_vector_byte_size + (ring_dimension * max_byte_size)
        let expected_total_size =
            4 + bit_vector_byte_size + (ring_dimension * max_byte_size as usize);
        assert_eq!(bytes.len(), expected_total_size, "Incorrect total byte size");

        // Check that the structure is as expected
        // Verify bit vector section exists
        let bit_vector = &bytes[4..4 + bit_vector_byte_size];
        assert_eq!(bit_vector.len(), bit_vector_byte_size, "Bit vector size is incorrect");

        // Verify coefficient values section exists
        let coeffs_section = &bytes[4 + bit_vector_byte_size..];
        assert_eq!(coeffs_section.len(), ring_dimension, "Coefficient section size is incorrect");

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
}
