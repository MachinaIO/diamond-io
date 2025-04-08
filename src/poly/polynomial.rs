use super::PolyElem;
use crate::poly::params::PolyParams;
use itertools::Itertools;
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
    + Send
    + Sync
{
    type Elem: PolyElem;
    type Params: PolyParams<Modulus = <Self::Elem as PolyElem>::Modulus>;
    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Self;
    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self;
    fn from_decomposed(params: &Self::Params, decomposed: &[Self]) -> Self;
    fn from_bytes(params: &Self::Params, bytes: &[u8]) -> Self {
        let log_q_bytes = params.modulus_bits().div_ceil(8);
        let dim = params.ring_dimension() as usize;
        debug_assert_eq!(bytes.len(), log_q_bytes * dim);
        let coeffs = bytes
            .chunks_exact(log_q_bytes)
            .map(|chunk| Self::Elem::from_bytes(&params.modulus(), chunk))
            .collect_vec();
        Self::from_coeffs(params, &coeffs)
    }
    fn from_compact_bytes(params: &Self::Params, bytes: &[u8]) -> Self;
    fn coeffs(&self) -> Vec<Self::Elem>;
    fn const_zero(params: &Self::Params) -> Self;
    fn const_one(params: &Self::Params) -> Self;
    fn const_minus_one(params: &Self::Params) -> Self;
    fn const_power_of_two(params: &Self::Params, k: usize) -> Self;
    fn const_rotate_poly(params: &Self::Params, shift: usize) -> Self {
        let zero = Self::const_zero(params);
        let mut coeffs = zero.coeffs();
        coeffs[shift] = Self::Elem::one(&params.modulus());
        Self::from_coeffs(params, &coeffs)
    }
    fn const_max(params: &Self::Params) -> Self;
    fn extract_highest_bits(&self) -> Vec<bool> {
        let mut bits = Vec::with_capacity(self.coeffs().len());
        for elem in self.coeffs() {
            bits.push(elem.extract_highest_bits());
        }
        bits
    }
    fn decompose_bits(&self, params: &Self::Params) -> Vec<Self>;
    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for elem in self.coeffs() {
            bytes.extend_from_slice(&elem.to_bytes());
        }
        bytes
    }
    fn to_compact_bytes(&self) -> Vec<u8>;
    fn to_bool_vec(&self) -> Vec<bool>;
    fn extract_bits_with_threshold(&self, params: &Self::Params) -> Vec<bool>;
}
