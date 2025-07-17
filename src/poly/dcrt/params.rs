use num_bigint::BigUint;
use num_traits::Num;
use openfhe::ffi;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, sync::Arc};

use crate::poly::PolyParams;

#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct DCRTPolyParams {
    /// polynomial ring dimension
    ring_dimension: u32,
    /// size of the tower
    crt_depth: usize,
    /// number of bits of each tower's modulus
    crt_bits: usize,
    /// ring modulus
    #[serde(skip)]
    modulus: Arc<BigUint>,
    /// bit size of the base for the gadget vector and decomposition
    base_bits: u32,
}

impl Debug for DCRTPolyParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyParams")
            .field("modulus", &self.modulus)
            .field("ring_dimension", &self.ring_dimension())
            .field("crt_depth", &self.crt_depth())
            .field("crt_bits", &self.crt_bits())
            .field("base_bits", &self.base_bits)
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

    fn base_bits(&self) -> u32 {
        self.base_bits
    }

    fn modulus_bits(&self) -> usize {
        self.modulus.bits() as usize
    }

    fn modulus_digits(&self) -> usize {
        self.crt_bits.div_ceil(self.base_bits as usize) * self.crt_depth
    }

    fn to_crt(&self) -> Vec<Arc<BigUint>> {
        ffi::GenCRTBasis(self.ring_dimension, self.crt_depth, self.crt_bits)
            .into_iter()
            .map(|m| Arc::new(BigUint::from_str_radix(&m, 10).expect("invalid CRT modulus")))
            .collect()
    }

    fn from_crt(crt: Vec<Self>) -> Self {
        todo!()
    }
}

impl Default for DCRTPolyParams {
    /// **note**  these parameters are insecure and only for test purpose
    fn default() -> Self {
        Self::new(4, 2, 17, 1)
    }
}

impl DCRTPolyParams {
    pub fn new(ring_dimension: u32, crt_depth: usize, crt_bits: usize, base_bits: u32) -> Self {
        // assert that ring_dimension is a power of 2
        assert!(ring_dimension.is_power_of_two(), "ring_dimension must be a power of 2");
        let modulus = ffi::GenModulus(ring_dimension, crt_depth, crt_bits);
        Self {
            ring_dimension,
            crt_depth,
            crt_bits,
            modulus: Arc::new(BigUint::from_str_radix(&modulus, 10).expect("invalid string")),
            base_bits,
        }
    }

    pub fn crt_depth(&self) -> usize {
        self.crt_depth
    }

    pub fn crt_bits(&self) -> usize {
        self.crt_bits
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_crt_decomposition() {
        let ring_dimension = 16;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;

        let params = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        let crt = params.to_crt();
        assert_eq!(crt.len(), crt_depth, "unexpected number of CRT primes");

        for (i, q) in crt.iter().enumerate() {
            assert_eq!(
                q.bits() as usize,
                crt_bits,
                "q[{}] has bit‑length {}, expected {}",
                i,
                q.bits(),
                crt_bits
            );
        }

        let product: BigUint = crt.iter().fold(BigUint::from(1u32), |acc, q| acc * &**q);
        assert_eq!(product, *params.modulus(), "product of q_i should equal Q");
    }

    #[test]
    fn test_params_initiation_ring_dimension() {
        let ring_dimension = 16;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), ring_dimension);
        assert_eq!(p.modulus_bits(), 204);
        assert_eq!(p.base_bits(), base_bits);

        let ring_dimension = 2;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), 2);
        assert_eq!(p.modulus_bits(), 204);

        let ring_dimension = 1;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), 1);
        assert_eq!(p.modulus_bits(), 204);
    }

    #[test]
    fn test_params_initiation_crt_depth() {
        let ring_dimension = 16;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), ring_dimension);
        assert_eq!(p.modulus_bits() as u32, (crt_depth * crt_bits) as u32);

        let ring_dimension = 16;
        let crt_depth = 5;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), ring_dimension);
        assert_eq!(p.modulus_bits() as u32, (crt_depth * crt_bits) as u32);

        let ring_dimension = 16;
        let crt_depth = 6;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), ring_dimension);
        assert_eq!(p.modulus_bits() as u32, (crt_depth * crt_bits) as u32);

        let ring_dimension = 16;
        let crt_depth = 7;
        let crt_bits = 20;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.ring_dimension(), ring_dimension);
        assert_eq!(p.modulus_bits() as u32, (crt_depth * crt_bits) as u32);
    }

    #[test]
    fn test_params_initiation_base() {
        let ring_dimension = 16;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.base_bits(), base_bits);

        let ring_dimension = 16;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 4;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.base_bits(), base_bits);

        let ring_dimension = 16;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 20;
        let p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits);
        assert_eq!(p.base_bits(), base_bits);
    }

    #[test]
    #[should_panic(expected = "ring_dimension must be a power of 2")]
    fn test_params_initiation_non_power_of_two() {
        let ring_dimension = 20;
        let crt_depth = 4;
        let crt_bits = 51;
        let base_bits = 1;
        let _p = DCRTPolyParams::new(ring_dimension, crt_depth, crt_bits, base_bits); // This should
                                                                                      // panic
    }
}
