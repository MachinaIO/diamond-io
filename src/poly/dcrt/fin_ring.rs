use std::{
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};

use crate::poly::PolyElem;
use num_bigint::{BigInt, BigUint};

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct FinRingElem {
    value: BigUint,
    modulus: Arc<BigUint>,
}

impl FinRingElem {
    pub fn new<V: Into<BigInt>>(value: V, modulus: Arc<BigUint>) -> Self {
        let value = value.into().to_biguint().unwrap();
        let reduced_value =
            if &value < modulus.as_ref() { value.clone() } else { value % modulus.as_ref() };
        Self { value: reduced_value, modulus }
    }

    pub fn value(&self) -> &BigUint {
        &self.value
    }

    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }

    pub fn modulus_switch(&self, new_modulus: Arc<BigUint>) -> Self {
        let value =
            ((&self.value * new_modulus.as_ref()) / self.modulus.as_ref()) % new_modulus.as_ref();
        Self { value, modulus: new_modulus }
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
        let max_minus_one = modulus.as_ref() - &BigUint::from(1u8);
        Self::new(max_minus_one, modulus.clone())
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
            return self;
        }
        Self::new(self.modulus.as_ref() - &self.value, self.modulus)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup() -> Arc<BigUint> {
        Arc::new(BigUint::from(17u8))
    }

    #[test]
    fn test_new_and_accessors() {
        let modulus = setup();

        // Test basic construction
        let elem = FinRingElem::new(5, modulus.clone());
        assert_eq!(elem.value(), &BigUint::from(5u8));
        assert_eq!(elem.modulus(), &BigUint::from(17u8));

        // Test reduction
        let elem = FinRingElem::new(20, modulus.clone()); // 20 ≡ 3 (mod 17)
        assert_eq!(elem.value(), &BigUint::from(3u8));
    }

    #[test]
    fn test_poly_elem_traits() {
        let modulus = setup();

        // Test zero, one, and minus_one
        let zero = FinRingElem::zero(&modulus);
        let one = FinRingElem::one(&modulus);
        let minus_one = FinRingElem::minus_one(&modulus);

        assert_eq!(zero.value(), &BigUint::from(0u8));
        assert_eq!(one.value(), &BigUint::from(1u8));
        assert_eq!(minus_one.value(), &BigUint::from(16u8)); // -1 ≡ 16 (mod 17)

        // Test extract_highest_bits
        let small = FinRingElem::new(3, modulus.clone());
        let large = FinRingElem::new(15, modulus.clone());
        assert!(small.extract_highest_bits()); // 3 < 17/2
        assert!(!large.extract_highest_bits()); // 15 > 17/2
    }

    #[test]
    fn test_arithmetic_operations() {
        let modulus = setup();

        // Addition tests
        let a = FinRingElem::new(5, modulus.clone());
        let b = FinRingElem::new(7, modulus.clone());

        assert_eq!((a.clone() + b.clone()).value(), &BigUint::from(12u8));
        assert_eq!((a.clone() + &b).value(), &BigUint::from(12u8));

        let mut c = a.clone();
        c += b.clone();
        assert_eq!(c.value(), &BigUint::from(12u8));

        let mut d = a.clone();
        d += &b;
        assert_eq!(d.value(), &BigUint::from(12u8));

        // Multiplication tests
        let prod = a.clone() * b.clone();
        assert_eq!(prod.value(), &BigUint::from(35u8 % 17u8));

        let mut e = a.clone();
        e *= b.clone();
        assert_eq!(e.value(), &BigUint::from(35u8 % 17u8));

        // Subtraction tests
        let diff = a.clone() - b.clone();
        assert_eq!(diff.value(), &BigUint::from(15u8)); // 5 - 7 ≡ -2 ≡ 15 (mod 17)

        let mut f = a.clone();
        f -= b.clone();
        assert_eq!(f.value(), &BigUint::from(15u8));

        // Negation test
        let neg_a = -a.clone();
        assert_eq!(neg_a.value(), &BigUint::from(12u8)); // -5 ≡ 12 (mod 17)

        // Test zero negation
        let zero = FinRingElem::zero(&modulus);
        let neg_zero = -zero;
        assert_eq!(neg_zero.value(), &BigUint::from(0u8));
    }

    #[test]
    fn test_arithmetic_with_reduction() {
        let modulus = setup();

        // Test addition with reduction
        let a = FinRingElem::new(15, modulus.clone());
        let b = FinRingElem::new(16, modulus.clone());
        assert_eq!((a + b).value(), &BigUint::from(14u8)); // 31 ≡ 14 (mod 17)

        // Test multiplication with reduction
        let c = FinRingElem::new(13, modulus.clone());
        let d = FinRingElem::new(15, modulus.clone());
        assert_eq!((c * d).value(), &BigUint::from(8u8)); // 195 ≡ 8 (mod 17)
    }

    #[test]
    fn test_arithmetic_modulus_switch() {
        let modulus = setup();

        let a = FinRingElem::new(15, modulus);
        let new_modulus = Arc::new(BigUint::from(5u8));
        assert_eq!(a.modulus_switch(new_modulus).value(), &BigUint::from(4u8));
    }
}
