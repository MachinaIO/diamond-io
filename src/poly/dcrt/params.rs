use crate::poly::params::PolyParams;
use num_bigint::BigUint;
use num_traits::Num;
use openfhe::ffi::{self};
use std::{fmt::Debug, sync::Arc};

#[derive(Clone)]
pub struct DCRTPolyParams {
    ring_dimension: u32,
    size: usize,
    k_res: usize,
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
        self.ring_dimension
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
        Self::new(16, 4, 51)
    }
}

impl DCRTPolyParams {
    pub fn new(n: u32, size: usize, k_res: usize) -> Self {
        // assert that n is a power of 2
        assert!(n.is_power_of_two(), "n must be a power of 2");
        let modulus = ffi::GenModulus(n, size, k_res);
        Self {
            ring_dimension: n,
            size,
            k_res,
            modulus: Arc::new(BigUint::from_str_radix(&modulus, 10).unwrap()),
        }
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn k_res(&self) -> usize {
        self.k_res
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

        let n = 2;
        let size = 4;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), 2);
        assert_eq!(p.modulus_bits(), 204);

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
        assert_eq!(p.modulus_bits() as u32, (size * k_res) as u32);

        let n = 16;
        let size = 5;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, (size * k_res) as u32);

        let n = 16;
        let size = 6;
        let k_res = 51;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, (size * k_res) as u32);

        let n = 16;
        let size = 7;
        let k_res = 20;
        let p = DCRTPolyParams::new(n, size, k_res);
        assert_eq!(p.ring_dimension(), n);
        assert_eq!(p.modulus_bits() as u32, (size * k_res) as u32);
    }

    #[test]
    #[should_panic(expected = "n must be a power of 2")]
    fn test_params_initiation_non_power_of_two() {
        let n = 20;
        let size = 4;
        let k_res = 51;
        let _p = DCRTPolyParams::new(n, size, k_res); // This should panic

        let n = 0;
        let size = 4;
        let k_res = 51;
        let _p = DCRTPolyParams::new(n, size, k_res); // This should panic
    }
}
