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
    fn scalar_mul(&self, params: &P::Params, scalar: &P) -> Self;
}

impl<P: Poly> Evaluable<P> for P {
    fn scalar_mul(&self, _: &P::Params, scalar: &P) -> Self {
        self.clone() * scalar
    }
}
