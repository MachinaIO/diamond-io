use super::Poly;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
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
    + Storable
    + Loadable
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
    fn concat_columns(&self, others: &[Self]) -> Self;
    /// (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[Self]) -> Self;
    /// (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[Self]) -> Self;
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
    fn modulus_switch(&self, new_params: &<Self::P as Poly>::Params) -> Self;
}

pub trait Storable {
    fn store(&self, path: &Path);
}

pub trait Loadable: Sized {
    fn load(path: &Path) -> Self;
}

pub struct PolyMatrixFSManager {
    base_dir: PathBuf,
}

impl PolyMatrixFSManager {
    /// Create a new PolyMatrixFSManager with the specified base directory
    /// The directory will be created if it doesn't exist
    pub fn new(base_dir: &Path) -> Self {
        fs::create_dir_all(base_dir).unwrap();
        Self { base_dir: base_dir.to_path_buf() }
    }

    /// Generate a path for a file
    pub fn gen_path(&self, name: &str) -> PathBuf {
        self.base_dir.join(name)
    }

    /// Clean up the directory
    pub fn cleanup(&self) -> std::io::Result<()> {
        fs::remove_dir_all(&self.base_dir)
    }

    /// Store a matrix with a given name and immediately drop it from memory
    pub fn store_matrix<M: PolyMatrix>(&self, name: &str, matrix: M) -> PathBuf {
        let path = self.gen_path(name);
        matrix.store(&path);
        path
    }

    /// Load a matrix with a given name
    pub fn load_matrix<M: PolyMatrix>(&self, name: &str) -> M {
        let path = self.gen_path(name);
        M::load(&path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams};
    use std::fs;

    fn setup_test_dir() -> PathBuf {
        let test_dir =
            std::env::temp_dir().join(format!("diamond-io-test-{}", rand::random::<u32>()));
        fs::create_dir_all(&test_dir).unwrap();
        test_dir
    }

    #[test]
    fn test_basic_store_load() {
        let test_dir = setup_test_dir();
        let fs_manager = PolyMatrixFSManager::new(&test_dir);

        // Store matrices
        let params = DCRTPolyParams::default();
        let sizes = [(1, 1), (2, 2), (3, 3)];
        for (idx, (nrow, ncol)) in sizes.iter().enumerate() {
            let poly_vec = vec![vec![DCRTPoly::const_one(&params); *ncol]; *nrow];
            let matrix = DCRTPolyMatrix::from_poly_vec(&params, poly_vec);
            fs_manager.store_matrix(&format!("test_matrix_{}", idx), matrix);
        }

        // Load matrices
        for i in 0..sizes.len() {
            let loaded = fs_manager.load_matrix::<DCRTPolyMatrix>(&format!("test_matrix_{}", i));
            let expected = DCRTPolyMatrix::from_poly_vec(
                &params,
                vec![vec![DCRTPoly::const_one(&params); sizes[i].1]; sizes[i].0],
            );
            assert_eq!(loaded, expected);
        }

        fs::remove_dir_all(test_dir.clone()).unwrap();

        assert!(!test_dir.exists());
    }

    #[test]
    #[should_panic]
    fn test_load_nonexistent() {
        let test_dir = setup_test_dir();
        let fs_manager = PolyMatrixFSManager::new(&test_dir);

        // Should panic when loading non-existent matrix
        let _matrix = fs_manager.load_matrix::<DCRTPolyMatrix>("nonexistent");

        fs::remove_dir_all(test_dir).unwrap();
    }
}
