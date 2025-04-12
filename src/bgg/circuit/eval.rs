use crate::poly::{Poly, PolyElem, PolyParams};
use rayon::prelude::*;
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
    fn from_digits(params: &Self::Params, one: &Self, digits: &[u32]) -> Self;
}

impl<P: Poly> Evaluable for P {
    type Params = P::Params;

    fn rotate(&self, params: &Self::Params, shift: usize) -> Self {
        let mut coeffs = self.coeffs();
        coeffs.rotate_right(shift);
        Self::from_coeffs(params, &coeffs)
    }

    fn from_digits(params: &Self::Params, _: &Self, digits: &[u32]) -> Self {
        let coeffs: Vec<P::Elem> = digits
            .par_iter()
            .map(|&digit| <P::Elem as PolyElem>::constant(&params.modulus(), digit as u64))
            .collect();
        Self::from_coeffs(params, &coeffs)
    }
}
