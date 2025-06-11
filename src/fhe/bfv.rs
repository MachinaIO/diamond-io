use num_bigint::BigUint;
use num_integer::Integer;
use std::ops::Add;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler, FinRingElem},
    sampler::{DistType, PolyUniformSampler},
    Poly, PolyParams,
};

pub struct Bfv {
    params: DCRTPolyParams,
}

#[derive(Debug)]
pub struct BfvCipher {
    c_1: DCRTPoly,
    c_2: DCRTPoly,
}

impl Bfv {
    pub fn keygen(params: DCRTPolyParams) -> (Self, DCRTPoly) {
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sk = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist).entry(0, 0);
        (Self { params }, sk)
    }

    pub fn encrypt_ske(
        &self,
        message: DCRTPoly,
        t: BigUint,
        sigma: f64,
        sk: DCRTPoly,
    ) -> BfvCipher {
        let value = self.params.modulus().div_ceil(&t);
        let delta_elem = FinRingElem::new(value, self.params.modulus());
        let delta_m = message.scalar_mul(&self.params, delta_elem);
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let a =
            sampler_uniform.sample_uniform(&self.params, 1, 1, DistType::FinRingDist).entry(0, 0);
        let e = sampler_uniform
            .sample_uniform(&self.params, 1, 1, DistType::GaussDist { sigma })
            .entry(0, 0);
        let c_1 = sk * &a + delta_m + e;
        let c_2 = -a;

        BfvCipher { c_1, c_2 }
    }
}

impl BfvCipher {
    pub fn decrypt(self, params: &DCRTPolyParams, sk: DCRTPoly) -> DCRTPoly {
        let ct = self.c_1 + self.c_2 * sk;
        ct
    }
}

impl Add for BfvCipher {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let c_1 = self.c_1 + rhs.c_1;
        let c_2 = self.c_2 + rhs.c_2;
        Self { c_1, c_2 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::create_random_poly;

    #[test]
    fn test_bfv_add_t_2_example() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        let (bfv, sk) = Bfv::keygen(params.clone());

        let m_a =
            create_random_poly(&params).modulus_switch(&params, BigUint::from(10 as u64).into());
        println!("m_a {:?}", m_a.coeffs());
        let enc_a = bfv.encrypt_ske(m_a, BigUint::from(10 as u64), 4.578, sk.clone());

        let m_b =
            create_random_poly(&params).modulus_switch(&params, BigUint::from(10 as u64).into());
        println!("m_b {:?}", m_b.coeffs());
        let enc_b = bfv.encrypt_ske(m_b, BigUint::from(10 as u64), 4.578, sk.clone());

        /* Homomorphic */
        let enc_3 = enc_a + enc_b;

        /* Decryption */
        // let raw_add = m_a + m_b;
        // println!("expected = {:?}", raw_add);
        let dec = enc_3.decrypt(&params, sk);
        println!("actual = {:?}", dec.coeffs());
    }
}
