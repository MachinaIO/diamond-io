use crate::poly::Poly;
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
    type Params;
    fn scalar_mul(&self, params: &Self::Params, scalar: &P) -> Self;
}

impl<P: Poly> Evaluable<P> for P {
    type Params = P::Params;
    fn scalar_mul(&self, _: &Self::Params, scalar: &P) -> Self {
        self.clone() * scalar
    }
}
