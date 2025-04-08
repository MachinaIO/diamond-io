#[cfg(feature = "disk")]
pub mod disk;
#[cfg(feature = "memory")]
pub mod memory;

#[cfg(feature = "disk")]
pub use disk::DCRTPolyMatrix;
#[cfg(feature = "memory")]
pub use memory::DCRTPolyMatrix;
