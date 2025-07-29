use super::{Poly, PolyParams};
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
    path::Path,
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

    fn to_compact_bytes(&self) -> Vec<u8>;
    fn from_compact_bytes(params: &<Self::P as Poly>::Params, bytes: &[u8]) -> Self;
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
    fn entry(&self, i: usize, j: usize) -> Self::P;
    fn set_entry(&mut self, i: usize, j: usize, elem: Self::P);
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
    fn zero(params: &<Self::P as Poly>::Params, nrow: usize, ncol: usize) -> Self;
    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self;
    fn transpose(&self) -> Self;
    /// (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[&Self]) -> Self;
    /// (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[&Self]) -> Self;
    /// (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[&Self]) -> Self;
    fn tensor(&self, other: &Self) -> Self;
    fn unit_column_vector(params: &<Self::P as Poly>::Params, size: usize, index: usize) -> Self {
        let mut vec = vec![Self::P::const_zero(params); size];
        vec[index] = Self::P::const_one(params);
        Self::from_poly_vec_column(params, vec)
    }
    fn unit_row_vector(params: &<Self::P as Poly>::Params, size: usize, index: usize) -> Self {
        let mut coeffs = vec![Self::P::const_zero(params); size];
        coeffs[index] = Self::P::const_one(params);
        Self::from_poly_vec_row(params, coeffs)
    }
    /// Constructs a gadget matrix Gₙ
    ///
    /// Gadget vector g = (b^0, b^1, ..., b^{log_b(q)-1}),
    /// where g ∈ Z_q^{log_b(q)} and b is the base defined in `params`.
    ///
    /// Gₙ = Iₙ ⊗ gᵀ
    ///
    /// * `params` - Parameters describing the modulus, the base, and other ring characteristics.
    /// * `size` - The size of the identity block (n), dictating the final matrix dimensions.
    ///
    /// A matrix of dimension n×(n·log_b(q)), in which each block row is a scaled identity
    /// under the ring modulus.
    fn gadget_matrix(params: &<Self::P as Poly>::Params, size: usize) -> Self;
    fn decompose(&self) -> Self;
    fn modulus_switch(
        &self,
        new_modulus: &<<Self::P as Poly>::Params as PolyParams>::Modulus,
    ) -> Self;
    /// Performs the operation S * (identity ⊗ other)
    fn mul_tensor_identity(&self, other: &Self, identity_size: usize) -> Self;
    /// Performs the operation S * (identity ⊗ G^-1(other)),
    /// where G^-1(other) is bit decomposition of other matrix
    fn mul_tensor_identity_decompose(&self, other: &Self, identity_size: usize) -> Self;
    /// j is column and return decomposed matrix of target column
    fn get_column_matrix_decompose(&self, j: usize) -> Self;
    /// Reads a matrix of given rows and cols with id from files under the given directory.
    fn read_from_files<P: AsRef<Path> + Send + Sync>(
        params: &<Self::P as Poly>::Params,
        nrow: usize,
        ncol: usize,
        dir_path: P,
        id: &str,
    ) -> Self;
    /// Extract block entries for parallel processing (used by storage service)
    fn block_entries(
        &self,
        rows: std::ops::Range<usize>,
        cols: std::ops::Range<usize>,
    ) -> Vec<Vec<Self::P>>;
}
