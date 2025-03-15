use super::{BggPublicKey, Evaluable};
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
        self.vector.concat_columns(&others.iter().map(|x| &x.vector).collect_vec()[..])
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
        let first_term = self.vector.clone() * decomposed_b.clone();
        let second_term = other.vector.clone() * self.plaintext.as_ref().unwrap();
        let new_vector = first_term + second_term;
        let new_plaintext = match (self.plaintext.as_ref(), other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a.clone() * b),
            _ => None,
        };

        let new_pubkey = BggPublicKey {
            matrix: self.pubkey.matrix.clone() * decomposed_b,
            reveal_plaintext: self.pubkey.reveal_plaintext & other.pubkey.reveal_plaintext,
        };
        Self { vector: new_vector, pubkey: new_pubkey, plaintext: new_plaintext }
    }
}

impl<M: PolyMatrix> Evaluable<M::P> for BggEncoding<M> {
    type Params = <M::P as Poly>::Params;
    fn scalar_mul(&self, params: &<M::P as Poly>::Params, scalar: &M::P) -> Self {
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
    use crate::{
        bgg::{
            circuit::PolyCircuit,
            sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        },
        poly::dcrt::{
            params::DCRTPolyParams,
            poly::DCRTPoly,
            sampler::{hash::DCRTPolyHashSampler, uniform::DCRTPolyUniformSampler},
        },
        utils::{create_bit_random_poly, create_random_poly},
    };
    use keccak_asm::Keccak256;
    use std::sync::Arc;

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
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, enc_one.clone(), &[enc1.clone(), enc2.clone()]);

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
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with a Sub operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, enc_one.clone(), &[enc1.clone(), enc2.clone()]);

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
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a simple circuit with a Mul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, enc_one.clone(), &[enc1.clone(), enc2.clone()]);

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
        let reveal_plaintexts = [true; 1];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
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
        let result = circuit.eval(&params, enc_one.clone(), &[enc.clone()]);

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
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
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
        let result =
            circuit.eval(&params, enc_one.clone(), &[enc1.clone(), enc2.clone(), enc3.clone()]);

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
        let reveal_plaintexts = [true; 4];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
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
        let result = circuit.eval(
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

    // #[test]
    // fn test_encoding_fhe_poly_bits_mul_by_poly_circuit() {
    //     // Create parameters for testing
    //     let params = DCRTPolyParams::new(4, 5, 7);

    //     // Create samplers
    //     let key: [u8; 32] = rand::random();
    //     let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
    //     let bgg_pubkey_sampler = BGGPublicKeySampler::new(hash_sampler);
    //     let uniform_sampler = Arc::new(DCRTPolyUniformSampler::new());

    //     // Generate random tag for sampling
    //     let tag: u64 = rand::random();
    //     let tag_bytes = tag.to_le_bytes();

    //     // Create random public keys
    //     let pubkeys =
    //         bgg_pubkey_sampler.sample(&params, &tag_bytes, (params.modulus_bits() * 2) + 2);

    //     // Create secret
    //     let secret = create_random_poly(&params);

    //     // Create plaintexts
    //     // encrypt a polynomial m using a RLWE secret key encryption
    //     // c0 = a*s + e + m (where m is the plaintext polynomial)
    //     // c1 = -a
    //     let m = uniform_sampler.sample_poly(&params, &DistType::BitDist);
    //     let s = uniform_sampler.sample_poly(&params, &DistType::BitDist);
    //     let e = uniform_sampler.sample_poly(&params, &DistType::GaussDist { sigma: 0.0 }); //
    // todo: set error     let a = uniform_sampler.sample_poly(&params, &DistType::FinRingDist);
    //     let c0 = -a.clone();
    //     let c1 = a * s.clone() + e + m.clone();

    //     // k is a polynomial from bit distribution
    //     let k = uniform_sampler.sample_poly(&params, &DistType::BitDist);

    //     let c0_bits = c0.decompose(&params);
    //     let c1_bits = c1.decompose(&params);

    //     // plaintexts is the concatenation of 1, c0_bits, c1_bits, k
    //     let plaintexts =
    //         [vec![DCRTPoly::const_one(&params)], c0_bits.clone(), c1_bits.clone(),
    // vec![k.clone()]]             .concat();

    //     assert_eq!(plaintexts.len(), (params.modulus_bits() * 2) + 2);

    //     // Create encoding sampler and encodings
    //     let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler,
    // 0.0);     let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts,
    // true);     let enc_one = encodings[0].clone();

    //     assert_eq!(encodings.len(), plaintexts.len());

    //     // Input: c0_bits[0], ..., c0_bits[modulus_bits - 1], c1_bits[0], ...,
    // c1_bits[modulus_bits - 1], k     // Output: c0_bits[0] * k, ..., c0_bits[modulus_bits -
    // 1] * k, c1_bits[0] * k, ..., c1_bits[modulus_bits - 1] * k     let mut circuit =
    // PolyCircuit::<DCRTPoly>::new();     let inputs = circuit.input((params.modulus_bits() *
    // 2) + 1);

    //     let k_id = inputs[inputs.len() - 1];
    //     let output_ids = inputs
    //         .iter()
    //         .take(inputs.len() - 1)
    //         .map(|&input_id| circuit.mul_gate(input_id, k_id))
    //         .collect();

    //     circuit.output(output_ids);

    //     // Evaluate the circuit
    //     let result = circuit.eval(&params, enc_one.clone(), &encodings[1..]);

    //     // Expected result: c0_bits_eval * k, c1_bits_eval * k
    //     let c0_bits_eval_bgg = result[..params.modulus_bits()].to_vec();
    //     let c1_bits_eval_bgg = result[params.modulus_bits()..].to_vec();

    //     let mut c0_bits_eval = Vec::with_capacity(params.modulus_bits());
    //     let mut c1_bits_eval = Vec::with_capacity(params.modulus_bits());

    //     for i in 0..params.modulus_bits() {
    //         c0_bits_eval.push(c0_bits_eval_bgg[i].plaintext.as_ref().unwrap().clone());
    //         c1_bits_eval.push(c1_bits_eval_bgg[i].plaintext.as_ref().unwrap().clone());
    //     }

    //     let c0_eval = DCRTPoly::from_decomposed(&params, &c0_bits_eval);
    //     let c1_eval = DCRTPoly::from_decomposed(&params, &c1_bits_eval);

    //     // Verify the result
    //     assert_eq!(result.len(), params.modulus_bits() * 2);
    //     assert_eq!(c0_eval, c0.clone() * k.clone());
    //     assert_eq!(c1_eval, c1.clone() * k.clone());

    //     // Decrypt the result
    //     let plaintext = c1_eval + c0_eval * s;
    //     assert_eq!(plaintext, m * k);
    // }

    #[test]
    fn test_encoding_register_and_call_sub_circuit() {
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
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![create_random_poly(&params), create_random_poly(&params)];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();

        // Create a sub-circuit that performs addition and multiplication
        let mut sub_circuit = PolyCircuit::<DCRTPoly>::new();
        let sub_inputs = sub_circuit.input(2);

        // Add operation: enc1 + enc2
        let add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);

        // Mul operation: enc1 * enc2
        let mul_gate = sub_circuit.mul_gate(sub_inputs[0], sub_inputs[1]);

        // Set the outputs of the sub-circuit
        sub_circuit.output(vec![add_gate, mul_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::<DCRTPoly>::new();
        let main_inputs = main_circuit.input(2);

        // Register the sub-circuit and get its ID
        let sub_circuit_id = main_circuit.register_sub_circuit(sub_circuit);

        // Call the sub-circuit with the main circuit's inputs
        let sub_outputs =
            main_circuit.call_sub_circuit(sub_circuit_id, &[main_inputs[0], main_inputs[1]]);

        // Verify we got two outputs from the sub-circuit
        assert_eq!(sub_outputs.len(), 2);

        // Use the sub-circuit outputs for further computation
        // For example, subtract the multiplication result from the addition result
        let final_gate = main_circuit.sub_gate(sub_outputs[0], sub_outputs[1]);

        // Set the output of the main circuit
        main_circuit.output(vec![final_gate]);

        // Evaluate the main circuit
        let result = main_circuit.eval(&params, enc_one, &[enc1.clone(), enc2.clone()]);

        // Expected result: (enc1 + enc2) - (enc1 * enc2)
        let expected = (enc1.clone() + enc2.clone()) - (enc1.clone() * enc2.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }

    #[test]
    fn test_encoding_nested_sub_circuits() {
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
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);

        // Create secret and plaintexts
        let secret = create_bit_random_poly(&params);
        let plaintexts = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create encoding sampler and encodings
        let bgg_encoding_sampler = BGGEncodingSampler::new(&params, &secret, uniform_sampler, 0.0);
        let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
        let enc_one = encodings[0].clone();
        let enc1 = encodings[1].clone();
        let enc2 = encodings[2].clone();
        let enc3 = encodings[3].clone();

        // Create a scalar
        let scalar = create_random_poly(&params);

        // Create the innermost sub-circuit that performs multiplication
        let mut inner_circuit = PolyCircuit::<DCRTPoly>::new();
        let inner_inputs = inner_circuit.input(2);
        let mul_gate = inner_circuit.mul_gate(inner_inputs[0], inner_inputs[1]);
        inner_circuit.output(vec![mul_gate]);

        // Create a middle sub-circuit that uses the inner sub-circuit
        let mut middle_circuit = PolyCircuit::<DCRTPoly>::new();
        let middle_inputs = middle_circuit.input(3);

        // Register the inner circuit
        let inner_circuit_id = middle_circuit.register_sub_circuit(inner_circuit);

        // Call the inner circuit with the first two inputs
        let inner_outputs = middle_circuit
            .call_sub_circuit(inner_circuit_id, &[middle_inputs[0], middle_inputs[1]]);

        // Add the result of the inner circuit with the third input
        let add_gate = middle_circuit.add_gate(inner_outputs[0], middle_inputs[2]);
        middle_circuit.output(vec![add_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::<DCRTPoly>::new();
        let main_inputs = main_circuit.input(3);

        // Register the middle circuit
        let middle_circuit_id = main_circuit.register_sub_circuit(middle_circuit);

        // Call the middle circuit with all inputs
        let middle_outputs = main_circuit
            .call_sub_circuit(middle_circuit_id, &[main_inputs[0], main_inputs[1], main_inputs[2]]);

        // Use the output for a scalar multiplication
        let scalar_mul_gate = main_circuit.scalar_mul_gate(middle_outputs[0], scalar.clone());

        // Set the output of the main circuit
        main_circuit.output(vec![scalar_mul_gate]);

        // Evaluate the main circuit
        let result =
            main_circuit.eval(&params, enc_one, &[enc1.clone(), enc2.clone(), enc3.clone()]);

        // Expected result: ((enc1 * enc2) + enc3) * scalar
        let expected = ((enc1.clone() * enc2.clone()) + enc3.clone()).scalar_mul(&params, &scalar);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].vector, expected.vector);
        assert_eq!(result[0].pubkey.matrix, expected.pubkey.matrix);
        assert_eq!(result[0].plaintext.as_ref().unwrap(), expected.plaintext.as_ref().unwrap());
    }
}
