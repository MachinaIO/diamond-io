use std::{
    ops::{Add, Mul},
    sync::Arc,
};

use crate::poly::{params::PolyElemParams, PolyElem};
use num_bigint::{BigInt, BigUint};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct FinRingParams {
    modulus: Arc<BigUint>,
}

impl PolyElemParams for FinRingParams {
    fn modulus(&self) -> BigUint {
        self.modulus.as_ref().clone()
    }
}

impl FinRingParams {
    pub fn new(modulus: Arc<BigUint>) -> Self {
        Self { modulus }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct FinRing {
    value: BigUint,
    params: FinRingParams,
}

impl FinRing {
    pub fn new<V: Into<BigInt>>(value: V, params: FinRingParams) -> Self {
        let value = value.into().to_biguint().unwrap();
        let modulus = params.modulus();
        let reduced_value = if value < modulus { value.clone() } else { value % modulus };
        Self { value: reduced_value, params }
    }

    pub fn value(&self) -> &BigUint {
        &self.value
    }

    pub fn modulus(&self) -> &BigUint {
        &self.params.modulus
    }
}

impl PolyElem for FinRing {
    type Params = FinRingParams;

    fn zero(params: &Self::Params) -> Self {
        Self::new(0, params.clone())
    }

    fn one(params: &Self::Params) -> Self {
        Self::new(1, params.clone())
    }

    fn minus_one(params: &Self::Params) -> Self {
        let max_minus_one = params.modulus() - &BigUint::from(1u8);
        Self::new(max_minus_one, params.clone())
    }

    fn extract_highest_bits(&self) -> bool {
        self.value < self.modulus() / 2u8
    }
}

// ====== Arithmetic ======

impl Add for FinRing {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value, self.params)
    }
}

impl Mul for FinRing {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.value * rhs.value, self.params)
    }
}
