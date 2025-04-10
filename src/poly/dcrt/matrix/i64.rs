use super::mmap_matrix::MmapMatrix;
use crate::poly::{MatrixElem, MatrixParams};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct I64MatrixParams;

impl MatrixParams for I64MatrixParams {
    fn entry_size(&self) -> usize {
        std::mem::size_of::<i64>()
    }
}

impl MatrixElem for i64 {
    type Params = I64MatrixParams;

    fn zero(_: &Self::Params) -> Self {
        0
    }

    fn one(_: &Self::Params) -> Self {
        1
    }

    fn from_bytes_to_elem(_: &Self::Params, bytes: &[u8]) -> Self {
        let arr: [u8; 8] = bytes.try_into().expect("slice length must be 8");
        i64::from_le_bytes(arr)
    }

    fn as_elem_to_bytes(&self) -> Vec<u8> {
        self.to_le_bytes().into()
    }
}

pub type I64Matrix = MmapMatrix<i64>;
