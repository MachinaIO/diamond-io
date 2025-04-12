use super::circuit::Evaluable;
use crate::poly::{Poly, PolyMatrix};
use rayon::prelude::*;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BggPublicKey<M: PolyMatrix> {
    pub matrix: M,
    pub reveal_plaintext: bool,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn new(matrix: M, reveal_plaintext: bool) -> Self {
        Self { matrix, reveal_plaintext }
    }

    pub fn concat_matrix(&self, others: &[Self]) -> M {
        self.matrix.concat_columns(&others.par_iter().map(|x| &x.matrix).collect::<Vec<_>>()[..])
    }
}

impl<M: PolyMatrix> Add for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix: self.matrix + &other.matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Sub for BggPublicKey<M> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self - &other
    }
}

impl<M: PolyMatrix> Sub<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn sub(self, other: &Self) -> Self {
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix: self.matrix - &other.matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Mul for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        let decomposed = other.matrix.decompose();
        let matrix = self.matrix.clone() * decomposed;
        let reveal_plaintext = self.reveal_plaintext & other.reveal_plaintext;
        Self { matrix, reveal_plaintext }
    }
}

impl<M: PolyMatrix> Evaluable<M::P> for BggPublicKey<M> {
    type Params = <M::P as Poly>::Params;
    fn rotate(&self, params: &Self::Params, shift: usize) -> Self {
        let rotate_poly = <M::P>::const_rotate_poly(params, shift);
        let matrix = self.matrix.clone() * rotate_poly;
        Self { matrix, reveal_plaintext: self.reveal_plaintext }
    }

    fn from_bits(params: &Self::Params, one: &Self, bits: &[bool]) -> Self {
        let const_poly =
            <M::P as Evaluable<M::P>>::from_bits(params, &<M::P>::const_one(params), bits);
        let matrix = one.matrix.clone() * const_poly;
        Self { matrix, reveal_plaintext: one.reveal_plaintext }
    }

    fn scalar_mul(&self, scalar: &M::P) -> Self {
        let matrix = self.matrix.clone() * scalar;
        let reveal_plaintext = self.reveal_plaintext;
        Self { matrix, reveal_plaintext }
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use crate::{
        bgg::{
            circuit::{Evaluable, PolyCircuit},
            sampler::BGGPublicKeySampler,
        },
        poly::dcrt::{params::DCRTPolyParams, DCRTPoly, DCRTPolyHashSampler},
        utils::create_random_poly,
    };
    use keccak_asm::Keccak256;
    use std::sync::Arc;

    #[test]
    fn test_pubkey_add() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2, 0).0;
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], &[]);

        // Expected result
        let expected = pk1.clone() + pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_sub() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with a Sub operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2, 0).0;
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], &[]);

        // Expected result
        let expected = pk1.clone() - pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with a Mul operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2, 0).0;
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], &[]);

        // Expected result
        let expected = pk1.clone() * pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_scalar_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public key
        let reveal_plaintexts = [true; 1];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();

        // Create scalar
        let scalar = create_random_poly(&params);

        // Create a simple circuit with a ScalarMul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, scalar_inputs) = circuit.input(1, 1);
        let scalar_mul_gate = circuit.scalar_mul_gate(inputs[0], scalar_inputs[0]);
        circuit.output(vec![scalar_mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone()], &[scalar.clone()]);

        // Expected result
        let expected = pk1.scalar_mul(&scalar);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_circuit_operations() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();

        // Create a circuit: ((pk1 + pk2)^2) - pk3
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(3, 0).0;

        // pk1 + pk2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (pk1 + pk2)^2
        let square_gate = circuit.mul_gate(add_gate, add_gate);

        // ((pk1 + pk2)^2) - pk3
        let sub_gate = circuit.sub_gate(square_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone(), pk3.clone()], &[]);

        // Expected result: ((pk1 + pk2)-2) - pk3
        let expected = ((pk1.clone() + pk2.clone()) * (pk1.clone() + pk2.clone())) - pk3.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_complex_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 4];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();
        let pk4 = pubkeys[4].clone();

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = pk1 + pk2, b = pk3 * pk4
        // Level 2: c = a * b, d = pk1 - pk3
        // Level 3: e = c + d
        // Level 4: f = e * e
        // Output: f
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(4, 0).0;

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // pk1 + pk2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // pk3 * pk4

        // Level 2
        let c = circuit.mul_gate(a, b); // (pk1 + pk2) * (pk3 * pk4)
        let d = circuit.sub_gate(inputs[0], inputs[2]); // pk1 - pk3

        // Level 3
        let e = circuit.add_gate(c, d); // ((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3)

        // Level 4
        let f = circuit.mul_gate(e, e); // (((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3))^2

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval(
            &params,
            &pk_one,
            &[pk1.clone(), pk2.clone(), pk3.clone(), pk4.clone()],
            &[],
        );

        // Expected result: (((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3))^2
        let sum1 = pk1.clone() + pk2.clone();
        let prod1 = pk3.clone() * pk4.clone();
        let prod2 = sum1.clone() * prod1;
        let diff = pk1.clone() - pk3.clone();
        let sum2 = prod2 + diff;
        let expected = sum2.clone() * sum2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_register_and_call_sub_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 2];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a sub-circuit that performs addition and multiplication
        let mut sub_circuit = PolyCircuit::new();
        let sub_inputs = sub_circuit.input(2, 0).0;

        // Add operation: pk1 + pk2
        let add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);

        // Mul operation: pk1 * pk2
        let mul_gate = sub_circuit.mul_gate(sub_inputs[0], sub_inputs[1]);

        // Set the outputs of the sub-circuit
        sub_circuit.output(vec![add_gate, mul_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::new();
        let main_inputs = main_circuit.input(2, 0).0;

        // Register the sub-circuit and get its ID
        let sub_circuit_id = main_circuit.register_sub_circuit(sub_circuit);

        // Call the sub-circuit with the main circuit's inputs
        let sub_outputs =
            main_circuit.call_sub_circuit(sub_circuit_id, &[main_inputs[0], main_inputs[1]], &[]);

        // Verify we got two outputs from the sub-circuit
        assert_eq!(sub_outputs.len(), 2);

        // Use the sub-circuit outputs for further computation
        // For example, subtract the multiplication result from the addition result
        let final_gate = main_circuit.sub_gate(sub_outputs[0], sub_outputs[1]);

        // Set the output of the main circuit
        main_circuit.output(vec![final_gate]);

        // Evaluate the main circuit
        let result = main_circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone()], &[]);

        // Expected result: (pk1 + pk2) - (pk1 * pk2)
        let expected = (pk1.clone() + pk2.clone()) - (pk1.clone() * pk2.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }

    #[test]
    fn test_pubkey_nested_sub_circuits() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let d = 3;
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler, d);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let reveal_plaintexts = [true; 3];
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();

        // Create the innermost sub-circuit that performs multiplication
        let mut inner_circuit = PolyCircuit::new();
        let inner_inputs = inner_circuit.input(2, 0).0;
        let mul_gate = inner_circuit.mul_gate(inner_inputs[0], inner_inputs[1]);
        inner_circuit.output(vec![mul_gate]);

        // Create a middle sub-circuit that uses the inner sub-circuit
        let mut middle_circuit = PolyCircuit::new();
        let middle_inputs = middle_circuit.input(3, 0).0;

        // Register the inner circuit
        let inner_circuit_id = middle_circuit.register_sub_circuit(inner_circuit);

        // Call the inner circuit with the first two inputs
        let inner_outputs = middle_circuit.call_sub_circuit(
            inner_circuit_id,
            &[middle_inputs[0], middle_inputs[1]],
            &[],
        );

        // Add the result of the inner circuit with the third input
        let add_gate = middle_circuit.add_gate(inner_outputs[0], middle_inputs[2]);
        middle_circuit.output(vec![add_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::new();
        let main_inputs = main_circuit.input(3, 0).0;

        // Register the middle circuit
        let middle_circuit_id = main_circuit.register_sub_circuit(middle_circuit);

        // Call the middle circuit with all inputs
        let middle_outputs = main_circuit.call_sub_circuit(
            middle_circuit_id,
            &[main_inputs[0], main_inputs[1], main_inputs[2]],
            &[],
        );

        // Use the output for square
        let square_gate = main_circuit.mul_gate(middle_outputs[0], middle_outputs[0]);

        // Set the output of the main circuit
        main_circuit.output(vec![square_gate]);

        // Evaluate the main circuit
        let result =
            main_circuit.eval(&params, &pk_one, &[pk1.clone(), pk2.clone(), pk3.clone()], &[]);

        // Expected result: ((pk1 * pk2) + pk3)^2
        let expected = ((pk1.clone() * pk2.clone()) + pk3.clone()) *
            ((pk1.clone() * pk2.clone()) + pk3.clone());

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
        assert_eq!(result[0].reveal_plaintext, expected.reveal_plaintext);
    }
}
