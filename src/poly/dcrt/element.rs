use crate::poly::PolyElem;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::Signed;
use std::{
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum FinRingElemError {
    #[error("Unknown plaintext for the left-hand input of multiplication: left-input index {0}, output index {1}")]
    UnknownPlaintextForMul(usize, usize),
}

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct FinRingElem {
    value: BigUint,
    modulus: Arc<BigUint>,
}

impl FinRingElem {
    pub fn new<V: Into<BigInt>>(value: V, modulus: Arc<BigUint>) -> Self {
        let value = value.into();
        let modulus_bigint = modulus.as_ref().to_bigint().unwrap();
        let value = if value.is_negative() {
            ((value % &modulus_bigint) + &modulus_bigint) % &modulus_bigint
        } else {
            value % &modulus_bigint
        };
        let reduced_value = value.to_biguint().unwrap();
        Self { value: reduced_value, modulus }
    }

    pub fn value(&self) -> &BigUint {
        &self.value
    }

    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }
}

impl PolyElem for FinRingElem {
    type Modulus = Arc<BigUint>;
    fn zero(modulus: &Self::Modulus) -> Self {
        Self::new(0, modulus.clone())
    }

    fn one(modulus: &Self::Modulus) -> Self {
        Self::new(1, modulus.clone())
    }

    fn minus_one(modulus: &Self::Modulus) -> Self {
        Self::new(modulus.as_ref() - &BigUint::from(1u8), modulus.clone())
    }

    fn extract_highest_bits(&self) -> bool {
        self.value < self.modulus() / 2u8
    }
}

// ====== Arithmetic ======

impl Add for FinRingElem {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a> Add<&'a FinRingElem> for FinRingElem {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self::new(self.value + &rhs.value, self.modulus)
    }
}

impl AddAssign for FinRingElem {
    fn add_assign(&mut self, rhs: Self) {
        *self = Self::new(&self.value + rhs.value, self.modulus.clone());
    }
}

impl<'a> AddAssign<&'a FinRingElem> for FinRingElem {
    fn add_assign(&mut self, rhs: &'a Self) {
        *self = Self::new(&self.value + &rhs.value, self.modulus.clone());
    }
}

impl Mul for FinRingElem {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a FinRingElem> for FinRingElem {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self::new(self.value * &rhs.value, self.modulus)
    }
}

impl MulAssign for FinRingElem {
    fn mul_assign(&mut self, rhs: Self) {
        *self = Self::new(&self.value * rhs.value, self.modulus.clone());
    }
}

impl<'a> MulAssign<&'a FinRingElem> for FinRingElem {
    fn mul_assign(&mut self, rhs: &'a Self) {
        *self = Self::new(&self.value * &rhs.value, self.modulus.clone());
    }
}

impl Sub for FinRingElem {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<'a> Sub<&'a FinRingElem> for FinRingElem {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self + (&-rhs.clone())
    }
}

impl SubAssign for FinRingElem {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl<'a> SubAssign<&'a FinRingElem> for FinRingElem {
    fn sub_assign(&mut self, rhs: &'a Self) {
        *self = self.clone() - rhs;
    }
}

impl Neg for FinRingElem {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.value == BigUint::from(0u8) {
            self
        } else {
            Self::new(self.modulus.as_ref() - &self.value, self.modulus)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_element_new() {
        let modulus = Arc::new(BigUint::from(17u8));
        let elem = FinRingElem::new(5, modulus.clone());
        assert_eq!(elem.value(), &BigUint::from(5u8));
        assert_eq!(elem.modulus(), &BigUint::from(17u8));
        let elem = FinRingElem::new(-5, modulus.clone());
        assert_eq!(elem.value(), &BigUint::from(12u8));
        assert_eq!(elem.modulus(), &BigUint::from(17u8));
        let elem = FinRingElem::new(-20, modulus.clone());
        assert_eq!(elem.value(), &BigUint::from(14u8));
        assert_eq!(elem.modulus(), &BigUint::from(17u8));
        let elem = FinRingElem::new(20, modulus.clone());
        assert_eq!(elem.value(), &BigUint::from(3u8));
    }

    #[test]
    fn test_element_zero() {
        let modulus = Arc::new(BigUint::from(17u8));
        let zero = FinRingElem::zero(&modulus);
        assert_eq!(zero.value(), &BigUint::from(0u8));
        assert_eq!(zero.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let zero = FinRingElem::zero(&modulus);
        assert_eq!(zero.value(), &BigUint::from(0u8));
        assert_eq!(zero.modulus(), modulus.as_ref());
    }

    #[test]
    fn test_element_one() {
        let modulus = Arc::new(BigUint::from(17u8));
        let one = FinRingElem::one(&modulus);
        assert_eq!(one.value(), &BigUint::from(1u8));
        assert_eq!(one.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let one = FinRingElem::one(&modulus);
        assert_eq!(one.value(), &BigUint::from(1u8));
        assert_eq!(one.modulus(), modulus.as_ref());
    }

    #[test]
    fn test_element_minus_one() {
        let modulus = Arc::new(BigUint::from(17u8));
        let minus_one = FinRingElem::minus_one(&modulus);
        assert_eq!(minus_one.value(), &BigUint::from(16u8));
        assert_eq!(minus_one.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let minus_one = FinRingElem::minus_one(&modulus);
        assert_eq!(minus_one.value(), &BigUint::from((10000 - 1) as usize));
        assert_eq!(minus_one.modulus(), modulus.as_ref());
    }

    #[test]
    fn test_element_higest_bit() {
        let modulus = Arc::new(BigUint::from(17u8));
        let small = FinRingElem::new(3, modulus.clone());
        let large = FinRingElem::new(15, modulus.clone());
        assert!(small.extract_highest_bits()); // 3 < 17/2
        assert!(!large.extract_highest_bits()); // 15 > 17/2

        let modulus = Arc::new(BigUint::from(10000usize));
        let small = FinRingElem::new(3, modulus.clone());
        let large = FinRingElem::new(15 + 10000, modulus.clone());
        assert!(small.extract_highest_bits()); // 3 < 10000/2
        assert!(large.extract_highest_bits()); // 15 + 10000 > 10000/2
    }

    #[test]
    fn test_element_add() {
        let modulus = Arc::new(BigUint::from(17u8));
        let a = FinRingElem::new(19, modulus.clone());
        let b = FinRingElem::new(16, modulus.clone());
        let c = a + b;
        assert_eq!(c.value(), &BigUint::from(1u8));
        assert_eq!(c.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let a = FinRingElem::new(19 + 10000, modulus.clone());
        let b = FinRingElem::new(16 + 10000, modulus.clone());
        let c = a + b;
        assert_eq!(c.value(), &BigUint::from(35u8));
        assert_eq!(c.modulus(), modulus.as_ref());
    }

    #[test]
    fn test_element_sub() {
        let modulus = Arc::new(BigUint::from(17u8));
        let a = FinRingElem::new(-1, modulus.clone());
        let b = FinRingElem::new(4, modulus.clone());
        let c = a - b;
        assert_eq!(c.value(), &BigUint::from(12u8));
        assert_eq!(c.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let a = FinRingElem::new(-19, modulus.clone());
        let b = FinRingElem::new(16 + 10000, modulus.clone());
        let c = a - b;
        assert_eq!(c.value(), &BigUint::from(9965usize));
        assert_eq!(c.modulus(), modulus.as_ref());
    }

    #[test]
    fn test_element_mul() {
        let modulus = Arc::new(BigUint::from(17u8));
        let a = FinRingElem::new(3, modulus.clone());
        let b = FinRingElem::new(5, modulus.clone());
        let c = a * b;
        assert_eq!(c.value(), &BigUint::from(15u8)); // 3 * 5 ≡ 15 (mod 17)
        assert_eq!(c.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let a = FinRingElem::new(200, modulus.clone());
        let b = FinRingElem::new(50, modulus.clone());
        let c = a * b;
        assert_eq!(c.value(), &BigUint::from(0usize)); // 200 * 50 ≡ 0 (mod 10000)
        assert_eq!(c.modulus(), modulus.as_ref());
    }

    #[test]
    fn test_element_neg() {
        let modulus = Arc::new(BigUint::from(17u8));
        let a = FinRingElem::new(5, modulus.clone());
        let neg_a = -a;
        assert_eq!(neg_a.value(), &BigUint::from(12u8)); // -5 ≡ 12 (mod 17)
        assert_eq!(neg_a.modulus(), modulus.as_ref());

        let modulus = Arc::new(BigUint::from(10000usize));
        let a = FinRingElem::new(200, modulus.clone());
        let neg_a = -a;
        assert_eq!(neg_a.value(), &BigUint::from(9800usize)); // -200 ≡ 9800 (mod 10000)
        assert_eq!(neg_a.modulus(), modulus.as_ref());
    }
}
