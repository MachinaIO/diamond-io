use itertools::Itertools;
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    path::Path,
};
use tokio;

use super::element::PolyElem;

pub trait PolyParams: Clone + Debug + PartialEq + Eq + Send + Sync {
    type Modulus: Debug + Clone;
    /// Returns the modulus value `q` used for polynomial coefficients in the ring `Z_q[x]/(x^n -
    /// 1)`.
    fn modulus(&self) -> Self::Modulus;
    /// A size of the base value used for a gadget vector and decomposition, i.e., `base =
    /// 2^base_bits`.
    fn base_bits(&self) -> u32;
    /// Fewest bits necessary to represent the modulus value `q`.
    fn modulus_bits(&self) -> usize;
    /// Fewest digits necessary to represent the modulus value `q` in the given base.
    fn modulus_digits(&self) -> usize;
    /// Returns the integer `n` that specifies the size of the polynomial ring used in this
    /// polynomial. Specifically, this is the degree parameter for the ring `Z_q[x]/(x^n - 1)`.
    fn ring_dimension(&self) -> u32;
}

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
    fn coeffs_digits(&self) -> Vec<u32> {
        self.coeffs()
            .iter()
            .map(|elem| {
                let u32s = elem.to_biguint().to_u32_digits();
                debug_assert!(u32s.len() < 2);
                if u32s.len() == 1 {
                    u32s[0]
                } else {
                    0
                }
            })
            .collect()
    }
    fn const_zero(params: &Self::Params) -> Self;
    fn const_one(params: &Self::Params) -> Self;
    fn const_minus_one(params: &Self::Params) -> Self;
    fn const_power_of_base(params: &Self::Params, k: usize) -> Self;
    fn const_int(params: &Self::Params, int: usize) -> Self;
    fn from_const_int_lsb(params: &Self::Params, int: usize) -> Self;
    fn const_rotate_poly(params: &Self::Params, shift: usize) -> Self {
        let zero = Self::const_zero(params);
        let mut coeffs = zero.coeffs();
        coeffs[shift] = Self::Elem::one(&params.modulus());
        Self::from_coeffs(params, &coeffs)
    }
    fn const_max(params: &Self::Params) -> Self;
    fn extract_bits_with_threshold(&self, params: &Self::Params) -> Vec<bool>;
    fn decompose_base(&self, params: &Self::Params) -> Vec<Self>;
    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        for elem in self.coeffs() {
            bytes.extend_from_slice(&elem.to_bytes());
        }
        bytes
    }
    fn to_bool_vec(&self) -> Vec<bool>;
    fn to_compact_bytes(&self) -> Vec<u8>;
    fn to_const_int(&self) -> usize;

    /// Reads a polynomial with id from files under the given directory.
    fn read_from_file<P: AsRef<Path> + Send + Sync>(
        params: &Self::Params,
        dir_path: P,
        id: &str,
    ) -> Self {
        let mut path = dir_path.as_ref().to_path_buf();
        path.push(format!("{id}.poly"));

        let bytes = std::fs::read(&path)
            .unwrap_or_else(|_| panic!("Failed to read polynomial file {path:?}"));

        Self::from_compact_bytes(params, &bytes)
    }

    /// Writes a polynomial with id to files under the given directory.
    fn write_to_file<P: AsRef<Path> + Send + Sync>(
        &self,
        dir_path: P,
        id: &str,
    ) -> impl std::future::Future<Output = ()> + Send {
        let mut path: std::path::PathBuf = dir_path.as_ref().to_path_buf();
        path.push(format!("{id}.poly"));

        let bytes = self.to_compact_bytes();
        async move {
            tokio::fs::write(&path, &bytes)
                .await
                .unwrap_or_else(|_| panic!("Failed to write polynomial file {path:?}"));
        }
    }
}
