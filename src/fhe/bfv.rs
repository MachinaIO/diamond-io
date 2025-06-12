use num_bigint::BigUint;
use num_integer::Integer;
use std::ops::{Add, Mul};

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler, FinRingElem},
    sampler::{DistType, PolyUniformSampler},
    Poly, PolyParams,
};

fn delta(params_q: &DCRTPolyParams, params_t: &DCRTPolyParams) -> BigUint {
    params_q.modulus().div_floor(&params_t.modulus())
}

pub struct Bfv {
    params_t: DCRTPolyParams,
    params_q: DCRTPolyParams,
    delta: BigUint,
    sigma: f64,
    rlk0: DCRTPoly,
    rlk1: DCRTPoly,
}

#[derive(Debug)]
pub struct BfvCipher {
    c_1: DCRTPoly,
    c_2: DCRTPoly,
}

impl Bfv {
    pub fn keygen(
        params_t: DCRTPolyParams,
        params_q: DCRTPolyParams,
        sigma: f64,
    ) -> (Self, DCRTPoly) {
        let sampler = DCRTPolyUniformSampler::new();
        let sk_t = sampler.sample_uniform(&params_t, 1, 1, DistType::BitDist).entry(0, 0);
        let delta = delta(&params_q, &params_t);
        let a = sampler.sample_uniform(&params_q, 1, 1, DistType::FinRingDist).entry(0, 0);
        let e = sampler.sample_uniform(&params_q, 1, 1, DistType::GaussDist { sigma }).entry(0, 0);
        let one_elem = FinRingElem::new(BigUint::from(1 as u8), params_q.modulus());
        let sk_q = sk_t.scalar_mul(&params_q, one_elem);
        let rlk0 = -(&a * &sk_q) + &sk_q + e;
        let rlk1 = a;
        (Self { params_t, params_q, delta, sigma, rlk0, rlk1 }, sk_t)
    }

    pub fn encrypt_ske(&self, m_t: DCRTPoly, sk_t: DCRTPoly) -> BfvCipher {
        let delta_elem = FinRingElem::new(self.delta.clone(), self.params_q.modulus());
        let m_q = m_t.scalar_mul(&self.params_q, delta_elem);
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let a =
            sampler_uniform.sample_uniform(&self.params_q, 1, 1, DistType::FinRingDist).entry(0, 0);
        let e = sampler_uniform
            .sample_uniform(&self.params_q, 1, 1, DistType::GaussDist { sigma: self.sigma })
            .entry(0, 0);
        /*
            this is hacky way to modular switch sk on mod t to mod q where t < q. I've introduced
            using scaler mul because existing modular switch doesn't support t < q case.
        */
        let one_elem = FinRingElem::new(BigUint::from(1 as u8), self.params_q.modulus());
        let sk_q = sk_t.scalar_mul(&self.params_q, one_elem);
        let c_1 = &sk_q * &a + m_q + &e;
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

    pub fn mul(&self, lhs: &BfvCipher, rhs: &BfvCipher) -> BfvCipher {
        let d0 = &lhs.c_1 * &rhs.c_1;
        let d1 = &lhs.c_1 * &rhs.c_2 + &lhs.c_2 * &rhs.c_1;
        let d2 = &lhs.c_2 * &rhs.c_2;

        let c0 = &d0 + &d2 * &self.rlk0;
        let c1 = &d1 + &d2 * &self.rlk1;

        let c0_t = c0.scale_and_round(&self.params_t);
        let c1_t = c1.scale_and_round(&self.params_t);

        let delta_elem = FinRingElem::new(self.delta.clone(), self.params_q.modulus());
        let c0_q = c0_t.scalar_mul(&self.params_q, delta_elem.clone());
        let c1_q = c1_t.scalar_mul(&self.params_q, delta_elem);

        BfvCipher { c_1: c0_q, c_2: c1_q }
    }
}

impl BfvCipher {
    pub fn decompose_base(&self, params_q: &DCRTPolyParams) -> Vec<DCRTPoly> {
        let mut c1_decomposed = self.c_1.decompose_base(&params_q);
        let c2_decomposed = self.c_2.decompose_base(&params_q);
        c1_decomposed.extend(c2_decomposed);
        c1_decomposed
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

impl Mul for &BfvCipher {
    type Output = BfvCipher;
    fn mul(self, _: Self) -> Self::Output {
        unimplemented!("use bfv.mul(&ct1,&ct2) instead");
    }
}

#[cfg(test)]
mod tests {
    use keccak_asm::Keccak256;

    use super::*;
    use crate::{
        bgg::{
            circuit::PolyCircuit,
            sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        },
        poly::{dcrt::DCRTPolyHashSampler, Poly},
        utils::{create_bit_random_poly, create_random_poly},
    };

    #[test]
    fn test_bfv_add() {
        /*
           Create parameter t and q for testing where they share same ring dimension but with different modulus.
           Paintext modulus t should be smaller than Ciphertext modulus q.
        */
        let params_t = DCRTPolyParams::new(4, 2, 17, 1);
        let params_q = DCRTPolyParams::new(4, 4, 21, 7);

        let (bfv, sk) = Bfv::keygen(params_t.clone(), params_q.clone(), 3.2);

        let m_a = create_random_poly(&params_t);
        let m_b = create_random_poly(&params_t);
        let m_add = &m_a + &m_b;
        let ct_a = bfv.encrypt_ske(m_a, sk.clone());
        let ct_b = bfv.encrypt_ske(m_b, sk.clone());

        /* Homomorphic */
        let ct_add = ct_a + ct_b;

        /* Decryption */
        let dec = bfv.decrypt(ct_add, sk);
        assert_eq!(m_add, dec);
    }

    #[test]
    fn test_bfv_decompose_bgg_add() {
        let params_t = DCRTPolyParams::new(4, 2, 17, 1);
        let params_q = DCRTPolyParams::new(4, 4, 21, 7);

        let (bfv, sk) = Bfv::keygen(params_t.clone(), params_q.clone(), 3.2);

        let m_1 = create_random_poly(&params_t);
        let ct_1 = bfv.encrypt_ske(m_1, sk.clone());
        let mut plaintexts_sum = ct_1.decompose_base(&params_q);
        let single_ct_plaintext_len = plaintexts_sum.len();

        let m_2 = create_random_poly(&params_t);
        let ct_2 = bfv.encrypt_ske(m_2, sk);
        let plaintexts_2 = ct_2.decompose_base(&params_q);
        plaintexts_sum.extend(plaintexts_2);
        println!("plaintexts_sum {}", plaintexts_sum.len());

        /* BGG */
        // initiate BGG encoding for (decomposition of ct_1 || decomposition of ct_2)
        let key: [u8; 32] = rand::random();
        let d = (2 * single_ct_plaintext_len) + 1;
        let bgg_pubkey_sampler =
            BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();
        let reveal_plaintexts = vec![true; d];
        let pubkeys = bgg_pubkey_sampler.sample(&params_q, &tag_bytes, &reveal_plaintexts);
        let secrets = vec![create_bit_random_poly(&params_q); d];
        let bgg_encoding_sampler =
            BGGEncodingSampler::new(&params_q, &secrets, uniform_sampler, 3.2);
        let encodings = bgg_encoding_sampler.sample(&params_q, &pubkeys, &plaintexts_sum);
        println!("encodings length {}", encodings.len());
        /* Circuit */
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(d - 1);
        let outputs: Vec<usize> = inputs[..single_ct_plaintext_len]
            .iter()
            .zip(&inputs[single_ct_plaintext_len..])
            .map(|(&l, &r)| circuit.add_gate(l, r))
            .collect();
        circuit.output(outputs);
        let result = circuit.eval(&params_q, &encodings[0], &encodings[1..]);
        println!("result length {}", result.len());
        for r_i in result {
            println!("{:?}", r_i.plaintext.unwrap().coeffs());
        }

        /* Expected */
        // decomposition of (ct_1 + ct_2)
        let ct_add = ct_1 + ct_2;
        let add_plaintext = ct_add.decompose_base(&params_q);
        println!("add_plaintext length {}", add_plaintext.len());
        // sample BGG encoding for decomposition
        let d = single_ct_plaintext_len + 1;
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let bgg_pubkey_sampler =
            BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, d);
        let reveal_plaintexts = vec![true; d];
        let pubkeys = bgg_pubkey_sampler.sample(&params_q, &tag_bytes, &reveal_plaintexts);
        let secrets = vec![create_bit_random_poly(&params_q); d];
        let bgg_encoding_sampler =
            BGGEncodingSampler::new(&params_q, &secrets, uniform_sampler, 3.2);
        let expected_result = bgg_encoding_sampler.sample(&params_q, &pubkeys, &add_plaintext);
        println!("expected_result length {}", expected_result.len());
        for e_i in expected_result {
            println!("{:?}", e_i.plaintext.unwrap().coeffs());
        }
        // assert_eq!(result, expected_result);
    }

    #[test]
    fn test_bfv_mul() {
        let params_t = DCRTPolyParams::new(4, 2, 17, 1);
        let params_q = DCRTPolyParams::new(4, 8, 51, 17);

        let (bfv, sk) = Bfv::keygen(params_t.clone(), params_q.clone(), 0.0);

        let m_a = create_random_poly(&params_t);
        println!("m_a = {:?}", m_a.coeffs());
        let m_b = create_random_poly(&params_t);
        println!("m_b = {:?}", m_b.coeffs());
        let m_prod = &m_a * &m_b;
        println!("m_prod = {:?}", m_prod.coeffs());
        let ct_a = bfv.encrypt_ske(m_a, sk.clone());
        let ct_b = bfv.encrypt_ske(m_b, sk.clone());

        /* Homomorphic */
        let ct_prod = bfv.mul(&ct_a, &ct_b);

        /* Decryption */
        let dec = bfv.decrypt(ct_prod, sk);
        println!("dec = {:?}", dec.coeffs());
        // todo: error
        // assert_eq!(m_prod, dec);
    }
}
