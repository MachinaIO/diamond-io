use super::{block_offsets, MmapMatrix, MmapMatrixElem, MmapMatrixParams};
use crate::{
    parallel_iter,
    poly::{
        dcrt::{DCRTPoly, DCRTPolyParams, FinRingElem},
        Poly, PolyMatrix, PolyParams,
    },
};
use itertools::Itertools;
use num_bigint::BigInt;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::ops::Range;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct I64MatrixParams;

impl MmapMatrixParams for I64MatrixParams {
    fn entry_size(&self) -> usize {
        std::mem::size_of::<i64>()
    }
}

impl MmapMatrixElem for i64 {
    type Params = I64MatrixParams;

    fn zero(_: &Self::Params) -> Self {
        0
    }

    fn one(_: &Self::Params) -> Self {
        1
    }

    fn from_bytes_to_elem(_: &Self::Params, bytes: &[u8]) -> Self {
        let mut arr = [0; 8];
        arr.copy_from_slice(bytes);
        i64::from_le_bytes(arr)
    }

    fn from_elem_to_bytes(&self) -> Vec<u8> {
        self.to_le_bytes().to_vec()
    }
}

pub type I64Matrix = MmapMatrix<i64>;
