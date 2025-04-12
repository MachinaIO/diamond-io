use crate::poly::{Poly, PolyElem, PolyParams};
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
};

pub trait Evaluable<P: Poly>:
    Debug
    + Clone
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
{
    type Params: Debug + Clone;
    fn rotate(&self, params: &Self::Params, shift: usize) -> Self;
    fn from_bits(params: &Self::Params, one: &Self, bits: &[bool]) -> Self;
    fn scalar_mul(&self, scalar: &P) -> Self;
}

impl<P: Poly> Evaluable<P> for P {
    type Params = P::Params;

    fn rotate(&self, params: &Self::Params, shift: usize) -> Self {
        let mut coeffs = self.coeffs();
        coeffs.rotate_right(shift);
        Self::from_coeffs(params, &coeffs)
    }

    fn from_bits(params: &Self::Params, _: &Self, bits: &[bool]) -> Self {
        let modulus = params.modulus();
        let one_elem = <P::Elem as PolyElem>::one(&modulus);
        let zero_elem = <P::Elem as PolyElem>::zero(&modulus);
        let coeffs: Vec<P::Elem> = bits
            .par_iter()
            .map(|&bit| if bit { one_elem.clone() } else { zero_elem.clone() })
            .collect();
        Self::from_coeffs(params, &coeffs)
    }

    fn scalar_mul(&self, scalar: &P) -> Self {
        self.clone() * scalar
    }
}
