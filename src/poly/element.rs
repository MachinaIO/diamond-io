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
    type Modulus: Debug + Clone;
    fn zero(modulus: &Self::Modulus) -> Self;
    fn one(modulus: &Self::Modulus) -> Self;
    fn minus_one(modulus: &Self::Modulus) -> Self;
    fn constant(modulus: &Self::Modulus, value: u64) -> Self;
    fn half_q(modulus: &Self::Modulus) -> Self;
    fn extract_highest_bits(&self) -> bool;
    fn modulus(&self) -> &Self::Modulus;
    fn to_biguint(&self) -> &num_bigint::BigUint;
}
