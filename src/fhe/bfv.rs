use num_bigint::BigUint;
use num_integer::Integer;
use std::ops::Add;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler, FinRingElem},
    sampler::{DistType, PolyUniformSampler},
    PolyParams,
};

fn delta(params_q: &DCRTPolyParams, params_t: &DCRTPolyParams) -> BigUint {
    params_q.modulus().div_floor(&params_t.modulus())
}

pub struct Bfv {
    params_t: DCRTPolyParams,
    params_q: DCRTPolyParams,
    delta: BigUint,
}

#[derive(Debug)]
pub struct BfvCipher {
    c_1: DCRTPoly,
    c_2: DCRTPoly,
}

impl Bfv {
    pub fn keygen(params_t: DCRTPolyParams, params_q: DCRTPolyParams) -> (Self, DCRTPoly) {
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sk = sampler_uniform.sample_uniform(&params_t, 1, 1, DistType::BitDist).entry(0, 0);
        let delta = delta(&params_q, &params_t);
        (Self { params_t, params_q, delta }, sk)
    }

    pub fn encrypt_ske(&self, m_t: DCRTPoly, sigma: f64, sk_t: DCRTPoly) -> BfvCipher {
        let delta_elem = FinRingElem::new(self.delta.clone(), self.params_q.modulus());
        let m_q = m_t.scalar_mul(&self.params_q, delta_elem);
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let a =
            sampler_uniform.sample_uniform(&self.params_q, 1, 1, DistType::FinRingDist).entry(0, 0);
        let e = sampler_uniform
            .sample_uniform(&self.params_q, 1, 1, DistType::GaussDist { sigma })
            .entry(0, 0);
        /*
            this is hacky way to modular switch sk on mod t to mod q where t < q. I've introduced
            using scaler mul because existing modular switch doesn't support t < q case.
        */
        let one_elem = FinRingElem::new(BigUint::from(1 as u8), self.params_q.modulus());
        let sk_q = sk_t.scalar_mul(&self.params_q, one_elem);
        let c_1 = sk_q * &a + m_q + e;
        let c_2 = -a;

        BfvCipher { c_1, c_2 }
    }

    pub fn decrypt(&self, ct: BfvCipher, sk_t: DCRTPoly) -> DCRTPoly {
        /*
            again, modular switch on sk from t to q
        */
        let one_elem = FinRingElem::new(BigUint::from(1 as u8), self.params_q.modulus());
        let sk_q = sk_t.scalar_mul(&self.params_q, one_elem);
        let ct = ct.c_1 + ct.c_2 * sk_q;
        ct.scale_and_round(&self.params_t)
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
    use crate::{poly::Poly, utils::create_random_poly};

    #[test]
    fn test_bfv_add_t_10_example() {
        /*
           Create parameter t and q for testing where they share same ring dimension but with different modulus.
           Paintext modulus t should be smaller than Ciphertext modulus q.
        */
        let params_t = DCRTPolyParams::new(4, 2, 17, 1);
        let params_q = DCRTPolyParams::new(4, 4, 21, 7);

        let (bfv, sk) = Bfv::keygen(params_t.clone(), params_q.clone());
        println!("sk {:?}", sk.coeffs());

        let m_a = create_random_poly(&params_t);
        let m_b = create_random_poly(&params_t);
        println!("m_a {:?}", m_a.coeffs());
        println!("m_b {:?}", m_b.coeffs());
        let raw_add = &m_a + &m_b;
        let enc_a = bfv.encrypt_ske(m_a, 0.0, sk.clone());
        let enc_b = bfv.encrypt_ske(m_b, 0.0, sk.clone());

        /* Homomorphic */
        let enc_3 = enc_a + enc_b;

        /* Decryption */
        println!("expected = {:?}", raw_add.coeffs());
        let dec = bfv.decrypt(enc_3, sk);
        println!("actual = {:?}", dec.coeffs());
        assert_eq!(raw_add, dec);
    }
}
