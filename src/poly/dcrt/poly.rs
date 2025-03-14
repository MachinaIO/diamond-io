use bytes::Bytes;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use super::{element::FinRingElem, params::DCRTPolyParams};
use crate::{
    parallel_iter,
    poly::{Poly, PolyElem, PolyParams},
};
use num_bigint::{BigInt, BigUint};
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPoly as DCRTPolyCxx},
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
            .iter()
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
        let coeffs = self.ptr_poly.GetCoefficients();
        let modulus = self.ptr_poly.GetModulus();
        parallel_iter!(coeffs)
            .map(|s| FinRingElem::from_str(&s, &modulus).expect("invalid string"))
            .collect()
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

    fn from_compact_bytes(params: &DCRTPolyParams, bytes: &Bytes, offset: usize) -> Self {
        let ring_dimension = params.ring_dimension() as usize;
        let modulus: BigUint = params.modulus().as_ref().clone();
        let modulus_big_int: BigInt =
            BigInt::from_bytes_le(num_bigint::Sign::Plus, &modulus.to_bytes_le());

        let byte_size = bytes.len() / ring_dimension;

        let mut coeffs = Vec::with_capacity(ring_dimension);

        for i in 0..ring_dimension {
            let start = i * byte_size;
            let end = start + byte_size;
            let value_bytes = &bytes[start..end];

            let value = (BigInt::from_bytes_le(num_bigint::Sign::Plus, value_bytes)
                - BigInt::from(offset))
                % modulus_big_int.clone();

            let coeff = FinRingElem::new(value, modulus.clone().into());

            coeffs.push(coeff);
        }

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
    fn decompose(&self, params: &Self::Params) -> Vec<Self> {
        let coeffs = self.coeffs();
        let bit_length = params.modulus_bits();
        parallel_iter!(0..bit_length)
            .map(|h| {
                let bit_coeffs: Vec<_> = coeffs
                    .iter()
                    .map(|j| {
                        let val = (j.value() >> h) & BigUint::from(1u32);
                        FinRingElem::new(val, params.modulus())
                    })
                    .collect();

                DCRTPoly::from_coeffs(params, &bit_coeffs)
            })
            .collect()
    }

    fn to_compact_bytes(&self, byte_size: usize, offset: usize) -> Bytes {
        let modulus = self.ptr_poly.GetModulus();
        let modulus_big: BigUint = BigUint::from_str(&modulus).unwrap();
        let modulus_bytes = modulus_big.to_bytes_le();

        assert!(
            byte_size <= modulus_bytes.len(),
            "byte_size must not be greater than modulus_bytes: {} > {}",
            byte_size,
            modulus_bytes.len()
        );

        let coeffs = self.coeffs();
        let mut bytes_data = Vec::new();
        for coeff in coeffs {
            let value = (coeff.value() + offset) % modulus_big.clone();
            let value_bytes = value.to_bytes_le();
            assert!(
                value_bytes.len() <= byte_size,
                "value_bytes exceeds the specified byte_size: {} > {}",
                value_bytes.len(),
                byte_size
            );

            if value_bytes.len() == byte_size {
                bytes_data.extend_from_slice(&value_bytes[0..byte_size]);
            } else {
                bytes_data.extend_from_slice(&value_bytes);
                bytes_data.extend(vec![0; byte_size - value_bytes.len()]);
            }
        }
        Bytes::from(bytes_data)
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

/// Implements $tr for all combinations of T and &T by delegating to the &T/&T implementation.
macro_rules! impl_binop_with_refs {
    ($T:ty => $tr:ident::$f:ident $($t:tt)*) => {
        impl $tr<$T> for $T {
            type Output = $T;

            #[inline]
            fn $f(self, rhs: $T) -> Self::Output {
                <&$T as $tr<&$T>>::$f(&self, &rhs)
            }
        }

        impl $tr<&$T> for $T {
            type Output = $T;

            #[inline]
            fn $f(self, rhs: &$T) -> Self::Output {
                <&$T as $tr<&$T>>::$f(&self, rhs)
            }
        }

        impl $tr<$T> for &$T {
            type Output = $T;

            #[inline]
            fn $f(self, rhs: $T) -> Self::Output {
                <&$T as $tr<&$T>>::$f(self, &rhs)
            }
        }

        impl $tr<&$T> for &$T {
            type Output = $T;

            #[inline]
            fn $f $($t)*
        }
    };
}

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
        let params = DCRTPolyParams::new(n, size, 51);
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
        let decomposed = poly.decompose(&params);
        assert_eq!(decomposed.len(), params.modulus_bits());
    }
}
