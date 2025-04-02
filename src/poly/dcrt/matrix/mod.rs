pub mod dcrt_poly_matrix;
pub mod i64_matrix;
pub mod memory_matrix;
pub mod mmap_matrix;

pub use dcrt_poly_matrix::DCRTPolyMatrix;
pub use i64_matrix::I64Matrix;
pub use mmap_matrix::{block_offsets, block_size, MmapMatrix, MmapMatrixElem, MmapMatrixParams};
