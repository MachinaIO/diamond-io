use crate::poly::params::PolyParams;
use num_bigint::BigUint;
use num_traits::Num;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParamsImpl},
};
use std::{fmt::Debug, sync::Arc};

#[derive(Clone)]
pub struct DCRTPolyParams {
    ptr_params: Arc<UniquePtr<ILDCRTParamsImpl>>,
    modulus: Arc<BigUint>,
}

impl Debug for DCRTPolyParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyParams")
            .field("modulus", &self.modulus())
            .field("ring_dimension", &self.ring_dimension())
            .finish()
    }
}

impl PolyParams for DCRTPolyParams {
    type Modulus = Arc<BigUint>;

    fn ring_dimension(&self) -> u32 {
        let ring_dimension = &self.ptr_params.as_ref().GetRingDimension();
        *ring_dimension
    }
    fn modulus(&self) -> Self::Modulus {
        self.modulus.clone()
    }
    fn modulus_bits(&self) -> usize {
        self.modulus().bits() as usize
    }
}

#[cfg(test)]
impl Default for DCRTPolyParams {
    fn default() -> Self {
        Self::new(4, 16, 51)
    }
}

impl DCRTPolyParams {
    pub fn new(n: u32, size: u32, k_res: u32) -> Self {
        let ptr_params = ffi::GenILDCRTParamsByOrderSizeBits(2 * n, size, k_res);
        let modulus = BigUint::from_str_radix(&ptr_params.GetModulus(), 10).unwrap();
        Self { ptr_params: ptr_params.into(), modulus: Arc::new(modulus) }
    }

    pub fn get_params(&self) -> &UniquePtr<ILDCRTParamsImpl> {
        &self.ptr_params
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_params_initiation() {
        let p = DCRTPolyParams::default();
        assert_eq!(p.ring_dimension(), 4);
        assert_eq!(p.modulus_bits(), 816);
    }
}
