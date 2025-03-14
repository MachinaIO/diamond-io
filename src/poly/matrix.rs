use super::{Poly, PolyParams};
use bytes::Bytes;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
};

pub trait PolyMatrix:
    Sized
    + Clone
    + Debug
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + Mul<Self::P, Output = Self>
    + for<'a> Mul<&'a Self::P, Output = Self>
    + Send
    + Sync
{
    type P: Poly;

    fn from_poly_vec(params: &<Self::P as Poly>::Params, vec: Vec<Vec<Self::P>>) -> Self;
    /// Creates a row vector (1 x n matrix) from a vector of n DCRTPoly elements.
    fn from_poly_vec_row(params: &<Self::P as Poly>::Params, vec: Vec<Self::P>) -> Self {
        // Wrap the vector in another vector to create a single row
        let wrapped_vec = vec![vec];
        Self::from_poly_vec(params, wrapped_vec)
    }
    /// Creates a column vector (n x 1 matrix) from a vector of DCRTPoly elements.
    fn from_poly_vec_column(params: &<Self::P as Poly>::Params, vec: Vec<Self::P>) -> Self {
        // Transform the vector into a vector of single-element vectors
        let wrapped_vec = vec.into_iter().map(|elem| vec![elem]).collect();
        Self::from_poly_vec(params, wrapped_vec)
    }
    fn from_compact_bytes(params: &<Self::P as Poly>::Params, bytes: Vec<Bytes>) -> Self;
    fn entry(&self, i: usize, j: usize) -> &Self::P;
    fn get_row(&self, i: usize) -> Vec<Self::P>;
    fn get_column(&self, j: usize) -> Vec<Self::P>;
    fn size(&self) -> (usize, usize);
    fn row_size(&self) -> usize {
        self.size().0
    }
    fn col_size(&self) -> usize {
        self.size().1
    }
    fn slice(
        &self,
        row_start: usize,
        row_end: usize,
        column_start: usize,
        column_end: usize,
    ) -> Self;
    fn slice_rows(&self, start: usize, end: usize) -> Self {
        let (_, columns) = self.size();
        self.slice(start, end, 0, columns)
    }
    fn slice_columns(&self, start: usize, end: usize) -> Self {
        let (rows, _) = self.size();
        self.slice(0, rows, start, end)
    }
    fn zero(params: &<Self::P as Poly>::Params, rows: usize, columns: usize) -> Self;
    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self;
    fn transpose(&self) -> Self;
    /// (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[&Self]) -> Self;
    /// (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[&Self]) -> Self;
    /// (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[&Self]) -> Self;
    fn tensor(&self, other: &Self) -> Self;
    /// Constructs a gadget matrix Gₙ
    ///
    /// Gadget vector g = (2^0, 2^1, ..., 2^{log(q)-1}),
    /// where g ∈ Z_q^{log(q)}.
    ///
    /// Gₙ = Iₙ ⊗ gᵀ
    ///
    /// * `params` - Parameters describing the modulus and other ring characteristics.
    /// * `size` - The size of the identity block (n), dictating the final matrix dimensions.
    ///
    /// A matrix of dimension n×(n·bit_length), in which each block row is a scaled identity
    /// under the ring modulus.
    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self;
    fn decompose(&self) -> Self;
    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self;
    /// Performs the operation S * (identity ⊗ other)
    fn mul_tensor_identity(&self, other: &Self, identity_size: usize) -> Self;
    fn to_compact_bytes(&self, byte_size: usize, offset: usize) -> Vec<Bytes>;
}
