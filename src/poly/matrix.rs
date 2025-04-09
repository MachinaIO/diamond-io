use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub trait MatrixParams: Debug + Clone + PartialEq + Eq + Send + Sync {
    fn entry_size(&self) -> usize;
}

pub trait MatrixElem:
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
    + Send
    + Sync
{
    type Params: MatrixParams;
    fn new_empty(params: &Self::Params) -> Self;
    fn zero(params: &Self::Params) -> Self;
    fn one(params: &Self::Params) -> Self;
    fn from_bytes_to_elem(params: &Self::Params, bytes: &[u8]) -> Self;
    fn as_elem_to_bytes(&self) -> Vec<u8>;
}
