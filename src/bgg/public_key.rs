use itertools::Itertools;
use std::ops::{Add, Mul, Sub};

use super::circuits::Evaluable;
use crate::poly::{Poly, PolyMatrix};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BggPublicKey<M: PolyMatrix> {
    pub matrix: M,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn new(matrix: M) -> Self {
        Self { matrix }
    }

    pub fn concat_matrix(&self, others: &[Self]) -> M {
        self.matrix.concat_columns(&others.iter().map(|x| x.matrix.clone()).collect_vec()[..])
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
        Self { matrix: self.matrix + &other.matrix }
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
        Self { matrix: self.matrix - &other.matrix }
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
        Self { matrix }
    }
}

impl<M: PolyMatrix> Evaluable<M::P> for BggPublicKey<M> {
    type Params = <M::P as Poly>::Params;
    fn scalar_mul(&self, params: &Self::Params, scalar: &M::P) -> Self {
        let gadget = M::gadget_matrix(params, 2);
        let scalared = gadget * scalar;
        let decomposed = scalared.decompose();
        let matrix = self.matrix.clone() * decomposed;
        Self { matrix }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgg::circuits::PolyCircuit;
    use crate::bgg::sampler::BGGPublicKeySampler;
    use crate::poly::dcrt::{
        params::DCRTPolyParams, poly::DCRTPoly, sampler::hash::DCRTPolyHashSampler,
        sampler::uniform::DCRTPolyUniformSampler,
    };
    use crate::poly::sampler::DistType;
    use keccak_asm::Keccak256;
    use std::sync::Arc;

    // Helper function to create a random polynomial using UniformSampler
    fn create_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
        let sampler = DCRTPolyUniformSampler::new();
        sampler.sample_poly(params, &DistType::FinRingDist)
    }

    #[test]
    fn test_pubkey_add() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, 3);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with an Add operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(&params, pk_one, &[pk1.clone(), pk2.clone()]);

        // Expected result
        let expected = pk1.clone() + pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
    }

    #[test]
    fn test_pubkey_sub() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, 3);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with a Sub operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(&params, pk_one, &[pk1.clone(), pk2.clone()]);

        // Expected result
        let expected = pk1.clone() - pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
    }

    #[test]
    fn test_pubkey_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, 3);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();

        // Create a simple circuit with a Mul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(&params, pk_one, &[pk1.clone(), pk2.clone()]);

        // Expected result
        let expected = pk1.clone() * pk2.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
    }

    #[test]
    fn test_pubkey_scalar_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public key
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, 2);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();

        // Create scalar
        let scalar = create_random_poly(&params);

        // Create a simple circuit with a ScalarMul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(1);
        let scalar_mul_gate = circuit.scalar_mul_gate(inputs[0], scalar.clone());
        circuit.output(vec![scalar_mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(&params, pk_one, &[pk1.clone()]);

        // Expected result
        let expected = pk1.scalar_mul(&params, &scalar);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
    }

    #[test]
    fn test_pubkey_circuit_operations() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, 4);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();

        // Create a scalar
        let scalar = create_random_poly(&params);

        // Create a circuit: ((pk1 + pk2) * scalar) - pk3
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(3);

        // pk1 + pk2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (pk1 + pk2) * scalar
        let scalar_mul_gate = circuit.scalar_mul_gate(add_gate, scalar.clone());

        // ((pk1 + pk2) * scalar) - pk3
        let sub_gate = circuit.sub_gate(scalar_mul_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval_poly_circuit(&params, pk_one, &[pk1.clone(), pk2.clone(), pk3.clone()]);

        // Expected result: ((pk1 + pk2) * scalar) - pk3
        let expected = ((pk1.clone() + pk2.clone()).scalar_mul(&params, &scalar)) - pk3.clone();

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
    }

    #[test]
    fn test_pubkey_complex_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a hash sampler and BGGPublicKeySampler to be reused
        let key: [u8; 32] = rand::random();
        let hash_sampler = Arc::new(DCRTPolyHashSampler::<Keccak256>::new(key));
        let bgg_sampler = BGGPublicKeySampler::new(hash_sampler);

        // Generate random tag for sampling
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        // Create random public keys
        let pubkeys = bgg_sampler.sample(&params, &tag_bytes, 5);
        let pk_one = pubkeys[0].clone();
        let pk1 = pubkeys[1].clone();
        let pk2 = pubkeys[2].clone();
        let pk3 = pubkeys[3].clone();
        let pk4 = pubkeys[4].clone();

        // Create a scalar
        let scalar = create_random_poly(&params);

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = pk1 + pk2, b = pk3 * pk4
        // Level 2: c = a * b, d = pk1 - pk3
        // Level 3: e = c + d
        // Level 4: f = e * scalar
        // Output: f
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(4);

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // pk1 + pk2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // pk3 * pk4

        // Level 2
        let c = circuit.mul_gate(a, b); // (pk1 + pk2) * (pk3 * pk4)
        let d = circuit.sub_gate(inputs[0], inputs[2]); // pk1 - pk3

        // Level 3
        let e = circuit.add_gate(c, d); // ((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3)

        // Level 4
        let f = circuit.scalar_mul_gate(e, scalar.clone()); // (((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3)) * scalar

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            pk_one,
            &[pk1.clone(), pk2.clone(), pk3.clone(), pk4.clone()],
        );

        // Expected result: (((pk1 + pk2) * (pk3 * pk4)) + (pk1 - pk3)) * scalar
        let sum1 = pk1.clone() + pk2.clone();
        let prod1 = pk3.clone() * pk4.clone();
        let prod2 = sum1.clone() * prod1;
        let diff = pk1.clone() - pk3.clone();
        let sum2 = prod2 + diff;
        let expected = sum2.scalar_mul(&params, &scalar);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].matrix, expected.matrix);
    }
}
