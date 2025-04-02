use crate::poly::{Poly, PolyElem, PolyParams};
use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
};

pub trait Evaluable:
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
}

impl<P: Poly> Evaluable for P {
    type Params = P::Params;

    fn rotate(&self, params: &Self::Params, shift: usize) -> Self {
        let mut coeffs = self.coeffs();
        coeffs.rotate_right(shift);
        Self::from_coeffs(params, &coeffs)
    }

    fn from_bits(params: &Self::Params, _: &Self, bits: &[bool]) -> Self {
        let poly = Self::const_zero(params);
        let one_elem = <P::Elem as PolyElem>::one(&params.modulus());
        let mut coeffs = poly.coeffs();
        for (coeff, &bit) in coeffs.iter_mut().zip(bits) {
            if bit {
                *coeff = one_elem.clone();
            }
        }
        Self::from_coeffs(params, &coeffs)
    }
}
