#[cfg(feature = "disk")]
pub mod disk;
#[cfg(not(feature = "disk"))]
pub mod memory;

#[cfg(feature = "disk")]
pub use disk::BaseMatrix;
#[cfg(not(feature = "disk"))]
pub use memory::BaseMatrix;
