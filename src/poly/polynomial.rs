use super::PolyElem;
use crate::poly::params::PolyParams;
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub trait Poly:
    Sized
    + Clone
    + Debug
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
{
    type Elem: PolyElem;
    type Params: PolyParams<Modulus = <Self::Elem as PolyElem>::Modulus>;
    type Error: Debug;
    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Result<Self, Self::Error>;
    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self;
    fn coeffs(&self) -> Result<Vec<Self::Elem>, Self::Error>;
    fn const_zero(params: &Self::Params) -> Self;
    fn const_one(params: &Self::Params) -> Self;
    fn const_minus_one(params: &Self::Params) -> Self;
    fn extract_highest_bits(&self) -> Result<Vec<bool>, Self::Error> {
        let mut bits = Vec::with_capacity(self.coeffs()?.len());
        for elem in self.coeffs()? {
            bits.push(elem.extract_highest_bits());
        }
        Ok(bits)
    }
}
