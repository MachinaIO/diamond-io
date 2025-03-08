use super::circuits::Evaluable;
use super::BggPublicKey;
use crate::poly::{Poly, PolyMatrix};
use itertools::Itertools;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone)]
pub struct BggEncoding<M: PolyMatrix> {
    pub vector: M,
    pub pubkey: BggPublicKey<M>,
    pub plaintext: Option<<M as PolyMatrix>::P>,
}

impl<M: PolyMatrix> BggEncoding<M> {
    pub fn new(
        vector: M,
        pubkey: BggPublicKey<M>,
        plaintext: Option<<M as PolyMatrix>::P>,
    ) -> Self {
        Self { vector, pubkey, plaintext }
    }

    pub fn concat_vector(&self, others: &[Self]) -> M {
        self.vector.concat_columns(&others.iter().map(|x| x.vector.clone()).collect_vec()[..])
    }
}

impl<M: PolyMatrix> Add for BggEncoding<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggEncoding<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        let vector = self.vector + &other.vector;
        let pubkey = self.pubkey + &other.pubkey;
        let plaintext = match (self.plaintext.as_ref(), other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a.clone() + b),
            _ => None,
        };
        Self { vector, pubkey, plaintext }
    }
}

impl<M: PolyMatrix> Sub for BggEncoding<M> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self - &other
    }
}

impl<M: PolyMatrix> Sub<&Self> for BggEncoding<M> {
    type Output = Self;
    fn sub(self, other: &Self) -> Self {
        let vector = self.vector - &other.vector;
        let pubkey = self.pubkey - &other.pubkey;
        let plaintext = match (self.plaintext.as_ref(), other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a.clone() - b),
            _ => None,
        };
        Self { vector, pubkey, plaintext }
    }
}

impl<M: PolyMatrix> Mul for BggEncoding<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggEncoding<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        if self.plaintext.is_none() {
            panic!("Unknown plaintext for the left-hand input of multiplication");
        }
        let decomposed_b = other.pubkey.matrix.decompose();
        let first_term = self.vector.clone() * decomposed_b;
        let second_term = other.vector.clone() * self.plaintext.as_ref().unwrap();
        let new_vector = first_term + second_term;
        let new_plaintext = match (self.plaintext.as_ref(), other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a.clone() * b),
            _ => None,
        };
        let new_pubkey = self.pubkey.clone() * &other.pubkey;
        Self { vector: new_vector, pubkey: new_pubkey, plaintext: new_plaintext }
    }
}

impl<M: PolyMatrix> Evaluable<M::P> for BggEncoding<M> {
    type Params = <M::P as Poly>::Params;

    fn scalar_mul(&self, params: &Self::Params, scalar: &M::P) -> Self {
        let gadget = M::gadget_matrix(params, 2);
        let scalared = gadget * scalar;
        let decomposed = scalared.decompose();
        let vector = self.vector.clone() * decomposed;
        let pubkey = self.pubkey.scalar_mul(params, scalar);
        let plaintext = self.plaintext.as_ref().map(|p| p.clone() * scalar);
        Self { vector, pubkey, plaintext }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgg::circuits::PolyCircuit;
    use crate::bgg::sampler::{BGGEncodingSampler, BGGPublicKeySampler};
    use crate::poly::dcrt::{
        params::DCRTPolyParams, poly::DCRTPoly, sampler::hash::DCRTPolyHashSampler,
        sampler::uniform::DCRTPolyUniformSampler,
    };
    use crate::poly::sampler::DistType;
    use crate::poly::PolyParams;
    use keccak_asm::Keccak256;
    use std::sync::Arc;

    // Helper function to create a random polynomial using UniformSampler
    fn create_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
        let sampler = DCRTPolyUniformSampler::new();
        sampler.sample_poly(params, &DistType::FinRingDist)
    }

    #[test]
    fn test_encoding_add() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, 3);

        // Create secret and plaintexts
        let secret = create_random_poly(&params);
        let plaintexts = vec![
            DCRTPoly::const_one(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval_poly_circuit(&params, enc_one.clone(), &[enc1.clone(), enc2.clone()]);

        // Expected result
        let expected = enc1.clone() + enc2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_sub() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, 3);

        // Create secret and plaintexts
        let secret = create_random_poly(&params);
        let plaintexts = vec![
            DCRTPoly::const_one(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with a Sub operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval_poly_circuit(&params, enc_one.clone(), &[enc1.clone(), enc2.clone()]);

        // Expected result
        let expected = enc1.clone() - enc2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, 3);

        // Create secret and plaintexts
        let secret = create_random_poly(&params);
        let plaintexts = vec![
            DCRTPoly::const_one(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with a Mul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval_poly_circuit(&params, enc_one.clone(), &[enc1.clone(), enc2.clone()]);

        // Expected result
        let expected = enc1.clone() * enc2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_scalar_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, 2);

        // Create secret and plaintexts
        let secret = create_random_poly(&params);
        let plaintexts = vec![DCRTPoly::const_one(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();
        let enc = encodings[1].clone();

        // Create scalar
        let scalar = create_random_poly(&params);

        // Create a simple circuit with a ScalarMul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(1);
        let scalar_mul_gate = circuit.scalar_mul_gate(inputs[0], scalar.clone());
        circuit.output(vec![scalar_mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(&params, enc_one.clone(), &[enc.clone()]);

        // Expected result
        let expected = enc.scalar_mul(&params, &scalar);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_circuit_operations() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, 4);

        // Create secret and plaintexts
        let secret = create_random_poly(&params);
        let plaintexts = vec![
            DCRTPoly::const_one(&params),
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();
        let enc3 = encodings[3].clone();

        // Create a scalar
        let scalar = create_random_poly(&params);

        // Create a circuit: ((enc1 + enc2) * scalar) - enc3
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(3);

        // enc1 + enc2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (enc1 + enc2) * scalar
        let scalar_mul_gate = circuit.scalar_mul_gate(add_gate, scalar.clone());

        // ((enc1 + enc2) * scalar) - enc3
        let sub_gate = circuit.sub_gate(scalar_mul_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            enc_one.clone(),
            &[enc1.clone(), enc2.clone(), enc3.clone()],
        );

        // Expected result: ((enc1 + enc2) * scalar) - enc3
        let expected = ((enc1.clone() + enc2.clone()).scalar_mul(&params, &scalar)) - enc3.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_complex_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, 5);

        // Create secret and plaintexts
        let secret = create_random_poly(&params);
        let plaintexts = vec![
            DCRTPoly::const_one(&params),
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();
        let enc3 = encodings[3].clone();
        let enc4 = encodings[4].clone();

        // Create a scalar
        let scalar = create_random_poly(&params);

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = enc1 + enc2, b = enc3 * enc4
        // Level 2: c = a * b, d = enc1 - enc3
        // Level 3: e = c + d
        // Level 4: f = e * scalar
        // Output: f
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(4);

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // enc1 + enc2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // enc3 * enc4

        // Level 2
        let c = circuit.mul_gate(a, b); // (enc1 + enc2) * (enc3 * enc4)
        let d = circuit.sub_gate(inputs[0], inputs[2]); // enc1 - enc3

        // Level 3
        let e = circuit.add_gate(c, d); // ((enc1 + enc2) * (enc3 * enc4)) + (enc1 - enc3)

        // Level 4
        let f = circuit.scalar_mul_gate(e, scalar.clone()); // (((enc1 + enc2) * (enc3 * enc4)) + (enc1 - enc3)) * scalar

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            enc_one.clone(),
            &[enc1.clone(), enc2.clone(), enc3.clone(), enc4.clone()],
        );

        // Expected result: (((enc1 + enc2) * (enc3 * enc4)) + (enc1 - enc3)) * scalar
        let sum1 = enc1.clone() + enc2.clone();
        let prod1 = enc3.clone() * enc4.clone();
        let prod2 = sum1.clone() * prod1;
        let diff = enc1.clone() - enc3.clone();
        let sum2 = prod2 + diff;
        let expected = sum2.scalar_mul(&params, &scalar);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_fhe_poly_bits_mul_by_poly_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::new(4, 5, 7);

        // Create samplers
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
        let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys =
            bgg_pubkey_sampler.sample(&params, &tag_bytes, (params.modulus_bits() * 2) + 2);

        // Create secret
        let secret = create_random_poly(&params);

        // Create plaintexts
        // encrypt a polynomial m using a RLWE secret key encryption
        // c0 = a*s + e + m (where m is the plaintext polynomial)
        // c1 = -a
        let m = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let s = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let e = uniform_sampler.sample_poly(&params, &DistType::GaussDist { sigma: 0.0 }); // todo: set error
        let a = uniform_sampler.sample_poly(&params, &DistType::FinRingDist);
        let c0 = -a.clone();
        let c1 = a * s.clone() + e + m.clone();

        // k is a polynomial from bit distribution
        let k = uniform_sampler.sample_poly(&params, &DistType::BitDist);

        let c0_bits = c0.decompose(&params);
        let c1_bits = c1.decompose(&params);

        // plaintexts is the concatenation of 1, c0_bits, c1_bits, k
        let plaintexts =
            [vec![DCRTPoly::const_one(&params)], c0_bits.clone(), c1_bits.clone(), vec![k.clone()]]
                .concat();

        assert_eq!(plaintexts.len(), (params.modulus_bits() * 2) + 2);

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts, true);
        let enc_one = encodings[0].clone();

        assert_eq!(encodings.len(), plaintexts.len());

        // Input: c0_bits[0], ..., c0_bits[modulus_bits - 1], c1_bits[0], ..., c1_bits[modulus_bits - 1], k
        // Output: c0_bits[0] * k, ..., c0_bits[modulus_bits - 1] * k, c1_bits[0] * k, ..., c1_bits[modulus_bits - 1] * k
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input((params.modulus_bits() * 2) + 1);

        let mut output_ids = Vec::with_capacity(c0_bits.len() + c1_bits.len());
        let k_id = inputs[inputs.len() - 1];
        for i in 0..inputs.len() - 1 {
            let output_id = circuit.mul_gate(inputs[i], k_id);
            output_ids.push(output_id);
        }

        circuit.output(output_ids);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(&params, enc_one.clone(), &encodings[1..]);

        // Expected result: c0_bits_eval * k, c1_bits_eval * k
        let c0_bits_eval_bgg = result[..params.modulus_bits()].to_vec();
        let c1_bits_eval_bgg = result[params.modulus_bits()..].to_vec();

        let mut c0_bits_eval = Vec::with_capacity(params.modulus_bits());
        let mut c1_bits_eval = Vec::with_capacity(params.modulus_bits());

        for i in 0..params.modulus_bits() {
            c0_bits_eval.push(c0_bits_eval_bgg[i].plaintext.as_ref().unwrap().clone());
            c1_bits_eval.push(c1_bits_eval_bgg[i].plaintext.as_ref().unwrap().clone());
        }

        let c0_eval = DCRTPoly::from_decomposed(&params, &c0_bits_eval);
        let c1_eval = DCRTPoly::from_decomposed(&params, &c1_bits_eval);

        // Verify the result
        assert_eq!(result.len(), params.modulus_bits() * 2);
        assert_eq!(c0_eval, c0.clone() * k.clone());
        assert_eq!(c1_eval, c1.clone() * k.clone());

        // Decrypt the result
        let plaintext = c1_eval + c0_eval * s;
        assert_eq!(plaintext, m * k);
    }
}
