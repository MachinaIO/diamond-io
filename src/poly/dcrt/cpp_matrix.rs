use openfhe::{
    cxx::UniquePtr,
    ffi::{GetMatrixCols, GetMatrixElement, GetMatrixRows, Matrix},
};

use super::DCRTPoly;

pub(crate) struct CppMatrix {
    pub(crate) inner: UniquePtr<Matrix>,
}

unsafe impl Send for CppMatrix {}
unsafe impl Sync for CppMatrix {}

impl CppMatrix {
    pub fn new(inner: UniquePtr<Matrix>) -> Self {
        CppMatrix { inner }
    }

    pub fn nrow(&self) -> usize {
        GetMatrixRows(&self.inner)
    }

    pub fn ncol(&self) -> usize {
        GetMatrixCols(&self.inner)
    }

    pub fn entry(&self, i: usize, j: usize) -> DCRTPoly {
        DCRTPoly::new(GetMatrixElement(&self.inner, i, j))
    }
}
