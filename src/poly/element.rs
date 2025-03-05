use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub trait PolyElem:
    Sized
    + Debug
    + Eq
    + Ord
    + Send
    + Sync
    + Clone
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
    type Modulus: Clone;
    fn zero(modulus: &Self::Modulus) -> Self;
    fn one(modulus: &Self::Modulus) -> Self;
    fn minus_one(modulus: &Self::Modulus) -> Self;
    fn extract_highest_bits(&self) -> bool;
}
