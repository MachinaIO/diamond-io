use super::{element::FinRingElem, params::DCRTPolyParams};
use crate::poly::{Poly, PolyParams};
use num_bigint::BigUint;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPoly as DCRTPolyCxx},
};
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};

#[derive(Clone, Debug)]
pub struct DCRTPoly {
    ptr_poly: Arc<UniquePtr<DCRTPolyCxx>>,
}

impl DCRTPoly {
    pub fn new(ptr_poly: UniquePtr<DCRTPolyCxx>) -> Self {
        Self { ptr_poly: ptr_poly.into() }
    }

    pub fn get_poly(&self) -> &UniquePtr<DCRTPolyCxx> {
        &self.ptr_poly
    }
}

impl Poly for DCRTPoly {
    type Elem = FinRingElem;
    type Params = DCRTPolyParams;

    fn coeffs(&self) -> Vec<Self::Elem> {
        let coeffs = self.ptr_poly.GetCoefficients();
        let mut result = Vec::with_capacity(coeffs.len());
        for s in coeffs.iter() {
            result.push(
                FinRingElem::from_str(s, &self.ptr_poly.GetModulus()).expect("invalid string"),
            );
        }
        result
    }

    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Self {
        let mut coeffs_cxx = Vec::with_capacity(coeffs.len());
        let modulus = params.modulus();
        for coeff in coeffs {
            let coeff_modulus = coeff.modulus();
            assert_eq!(coeff_modulus, modulus.as_ref());
            coeffs_cxx.push(coeff.value().to_string());
        }
        DCRTPoly::new(ffi::DCRTPolyGenFromVec(
            params.ring_dimension(),
            params.size(),
            params.k_res(),
            &coeffs_cxx,
        ))
    }

    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.ring_dimension(),
            params.size(),
            params.k_res(),
            &constant.value().to_string(),
        ))
    }

    fn const_zero(params: &Self::Params) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.ring_dimension(),
            params.size(),
            params.k_res(),
            &BigUint::ZERO.to_string(),
        ))
    }

    fn const_one(params: &Self::Params) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.ring_dimension(),
            params.size(),
            params.k_res(),
            &BigUint::from(1u32).to_string(),
        ))
    }

    fn const_minus_one(params: &Self::Params) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.ring_dimension(),
            params.size(),
            params.k_res(),
            &(params.modulus().as_ref() - BigUint::from(1u32)).to_string(),
        ))
    }
}

impl Add for DCRTPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a> Add<&'a DCRTPoly> for DCRTPoly {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        DCRTPoly::new(ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly))
    }
}

impl Mul for DCRTPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a DCRTPoly> for DCRTPoly {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        DCRTPoly::new(ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly))
    }
}

impl Sub for DCRTPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<'a> Sub<&'a DCRTPoly> for DCRTPoly {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self + rhs.clone().neg()
    }
}

impl Neg for DCRTPoly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        DCRTPoly::new(self.ptr_poly.Negate())
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

impl AddAssign for DCRTPoly {
    fn add_assign(&mut self, rhs: Self) {
        self.ptr_poly = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly).into();
    }
}

impl<'a> AddAssign<&'a DCRTPoly> for DCRTPoly {
    fn add_assign(&mut self, rhs: &'a Self) {
        self.ptr_poly = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly).into();
    }
}

impl MulAssign for DCRTPoly {
    fn mul_assign(&mut self, rhs: Self) {
        self.ptr_poly = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly).into();
    }
}

impl<'a> MulAssign<&'a DCRTPoly> for DCRTPoly {
    fn mul_assign(&mut self, rhs: &'a Self) {
        self.ptr_poly = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly).into();
    }
}

impl SubAssign for DCRTPoly {
    fn sub_assign(&mut self, rhs: Self) {
        let neg_rhs = rhs.neg();
        self.ptr_poly = ffi::DCRTPolyAdd(&neg_rhs.ptr_poly, &self.ptr_poly).into();
    }
}

impl<'a> SubAssign<&'a DCRTPoly> for DCRTPoly {
    fn sub_assign(&mut self, rhs: &'a Self) {
        let neg_rhs = rhs.clone().neg();
        self.ptr_poly = ffi::DCRTPolyAdd(&neg_rhs.ptr_poly, &self.ptr_poly).into();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::PolyParams;
    use rand::prelude::*;

    #[test]
    fn test_dcrtpoly_coeffs() {
        let mut rng = rand::rng();
        // todo: if x=0, n=1: libc++abi: terminating due to uncaught exception of type lbcrypto::OpenFHEException: /Users/piapark/Documents/GitHub/openfhe-development/src/core/include/math/nbtheory.h:l.156:ReverseBits(): msbb value not handled:0
        // todo: if x=1, n=2: value mismatch from_coeffs & coeffs
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
        let product = poly1.clone() * poly2.clone();

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
}
