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
            .field("modulus", &self.modulus)
            .field("ring_dimension", &self.ring_dimension())
            .finish()
    }
}

impl PolyParams for DCRTPolyParams {
    type Modulus = Arc<BigUint>;

    fn ring_dimension(&self) -> u32 {
        self.ptr_params.as_ref().GetRingDimension()
    }
    fn modulus(&self) -> Self::Modulus {
        self.modulus.clone()
    }
    fn modulus_bits(&self) -> usize {
        self.modulus.bits() as usize
    }
}

#[cfg(test)]
impl Default for DCRTPolyParams {
    fn default() -> Self {
        Self::new(4, 16, 51)
    }
}

impl DCRTPolyParams {
    /// 2 * n = order, size = depth, bits= k
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
    fn test_params_initiation_n() {
        let n = 16;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits(), 204);

        let n = 20;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        // ring dimension returning closest 2^n form
        assert_eq!(p.ring_dimension(), 16);
        assert_eq!(p.modulus_bits(), 204);

        let n = 2;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), 2);
        assert_eq!(p.modulus_bits(), 204);

        let n = 0;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), 0);
        assert_eq!(p.modulus_bits(), 0);

        let n = 1;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), 1);
        assert_eq!(p.modulus_bits(), 204);
    }

    #[test]
    fn test_params_initiation_size() {
        let n = 16;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, size * k_res);

        let n = 16;
        let size = 5;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, size * k_res);

        let n = 16;
        let size = 6;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, size * k_res);

        let n = 16;
        let size = 7;
        let k_res = 20;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, size * k_res);
    }
}
