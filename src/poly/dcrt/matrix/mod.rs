#[cfg(feature = "disk")]
pub mod disk;
#[cfg(feature = "memory")]
pub mod memory;

#[cfg(feature = "disk")]
pub use disk::{DCRTPolyMatrix, I64Matrix, I64MatrixParams};
#[cfg(feature = "memory")]
pub use memory::{DCRTPolyMatrix, I64Matrix, I64MatrixParams};
