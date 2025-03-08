use itertools::Itertools;

use super::*;
use crate::poly::Poly;
use std::collections::BTreeMap;

#[derive(Debug, Clone)]
pub struct PolyCircuit<P: Poly> {
    pub gates: BTreeMap<usize, PolyGate<P>>,
    pub output_ids: Vec<usize>,
    pub num_input: usize,
}

impl<P: Poly> Default for PolyCircuit<P> {
    fn default() -> Self {
        Self::new()
    }
}

impl<P: Poly> PolyCircuit<P> {
    pub fn new() -> Self {
        Self { gates: BTreeMap::new(), output_ids: vec![], num_input: 0 }
    }

    pub fn num_input(&self) -> usize {
        self.num_input
    }

    pub fn num_output(&self) -> usize {
        self.output_ids.len()
    }

    pub fn num_gates(&self) -> usize {
        self.gates.len()
    }

    pub fn input(&mut self, num_input: usize) -> Vec<usize> {
        #[cfg(debug_assertions)]
        assert_eq!(self.num_input, 0);
        self.gates.insert(0, PolyGate::new(0, PolyGateType::Input, vec![]));
        let mut input_gates = Vec::with_capacity(num_input);
        for idx in 1..(num_input + 1) {
            self.gates.insert(idx, PolyGate::new(idx, PolyGateType::Input, vec![]));
            input_gates.push(idx);
        }
        self.num_input += num_input;
        input_gates
    }

    pub fn output(&mut self, output_gates: Vec<usize>) {
        #[cfg(debug_assertions)]
        assert_eq!(self.output_ids.len(), 0);

        for gate_id in output_gates.iter() {
            self.output_ids.push(*gate_id);
        }
    }

    pub fn and_gate(&mut self, left: usize, right: usize) -> usize {
        self.mul_gate(left, right)
    }

    /// Computes the NOT gate using arithmetic inversion: `1 - x`.
    /// This operation assumes that `x` is restricted to binary values (0 or 1),
    /// meaning it should only be used with polynomials sampled from a bit distribution.
    /// The computation is achieved by subtracting `x` from 1 (i.e., `0 - x + 1`).
    pub fn not_gate(&mut self, input: usize) -> usize {
        self.sub_gate(0, input)
    }

    pub fn or_gate(&mut self, left: usize, right: usize) -> usize {
        let add = self.add_gate(left, right);
        let mul = self.mul_gate(left, right);
        self.sub_gate(add, mul) // A + B - A*B
    }

    /// Computes the NAND gate as `NOT(AND(left, right))`.
    /// This operation follows the same restriction as the NOT gate:
    /// `left` and `right` must be bit distribution (0 or 1)
    pub fn nand_gate(&mut self, left: usize, right: usize) -> usize {
        let and_result = self.and_gate(left, right);
        self.not_gate(and_result) // NOT AND
    }

    /// Computes the NOR gate as `NOT(OR(left, right))`.
    /// This operation follows the same restriction as the NOT gate:
    /// `left` and `right` must be bit distribution (0 or 1)
    pub fn nor_gate(&mut self, left: usize, right: usize) -> usize {
        let or_result = self.or_gate(left, right);
        self.not_gate(or_result) // NOT OR
    }

    pub fn xor_gate(&mut self, left: usize, right: usize) -> usize {
        let two = self.add_gate(0, 0);
        let mul = self.mul_gate(left, right);
        let two_mul = self.mul_gate(two, mul);
        let add = self.add_gate(left, right);
        self.sub_gate(add, two_mul) // A + B - 2*A*B
    }

    /// Computes the XNOR gate as `NOT(XOR(left, right))`.
    /// This operation follows the same restriction as the NOT gate:
    /// `left` and `right` must be bit distribution (0 or 1)
    pub fn xnor_gate(&mut self, left: usize, right: usize) -> usize {
        let xor_result = self.xor_gate(left, right);
        self.not_gate(xor_result) // NOT XOR
    }

    pub(crate) fn add_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Add)
    }

    pub(crate) fn sub_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Sub)
    }

    pub fn scalar_mul_gate(&mut self, input: usize, scalar: P) -> usize {
        self.new_gate_generic(vec![input], PolyGateType::ScalarMul(scalar))
    }

    pub(crate) fn mul_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Mul)
    }

    fn new_gate_generic(&mut self, input_gates: Vec<usize>, gate_type: PolyGateType<P>) -> usize {
        #[cfg(debug_assertions)]
        {
            assert_ne!(self.num_input, 0);
            assert_eq!(self.output_ids.len(), 0);
            assert_eq!(input_gates.len(), gate_type.num_input());
            for gate_id in input_gates.iter() {
                assert!(self.gates.contains_key(gate_id));
            }
        }
        let gate_id = self.gates.len();
        self.gates.insert(gate_id, PolyGate::new(gate_id, gate_type, input_gates));
        gate_id
    }

    pub fn eval_poly_circuit<E: Evaluable<P>>(
        self,
        params: &E::Params,
        one: E,
        input: &[E],
    ) -> Vec<E> {
        #[cfg(debug_assertions)]
        {
            assert_eq!(self.num_input(), input.len());
        }

        let mut wires = BTreeMap::new();
        wires.insert(0, one);
        for (idx, input) in input.iter().enumerate() {
            wires.insert(idx + 1, input.clone());
        }
        let output_gates = &self.output_ids.iter().map(|id| self.gates[id].clone()).collect_vec();
        for gate in output_gates.iter() {
            self.eval_poly_gate(params, &mut wires, gate);
        }
        let output = self.output_ids.iter().map(|id| wires.get(id).unwrap().clone()).collect_vec();
        output
    }

    fn eval_poly_gate<E: Evaluable<P>>(
        &self,
        params: &E::Params,
        wires: &mut BTreeMap<usize, E>,
        gate: &PolyGate<P>,
    ) {
        let input_ids = &gate.input_gates;
        for input_id in input_ids.iter() {
            if !wires.contains_key(input_id) {
                let input_gate = &self.gates[input_id];
                self.eval_poly_gate(params, wires, input_gate);
            }
        }
        match &gate.gate_type {
            PolyGateType::Input => {
                panic!("The wire for the input gate should be already inserted");
            }
            PolyGateType::Add => {
                let left = wires.get(&input_ids[0]).unwrap();
                let right = wires.get(&input_ids[1]).unwrap();
                let output = left.clone() + right;
                wires.insert(gate.gate_id, output);
            }
            PolyGateType::Sub => {
                let left = wires.get(&input_ids[0]).unwrap();
                let right = wires.get(&input_ids[1]).unwrap();
                let output = left.clone() - right;
                wires.insert(gate.gate_id, output);
            }
            PolyGateType::ScalarMul(scalar) => {
                let input = wires.get(&input_ids[0]).unwrap();
                let output = input.scalar_mul(params, scalar);
                wires.insert(gate.gate_id, output);
            }
            PolyGateType::Mul => {
                let left = wires.get(&input_ids[0]).unwrap();
                let right = wires.get(&input_ids[1]).unwrap();
                let output = left.clone() * right.clone();
                wires.insert(gate.gate_id, output);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{
            params::DCRTPolyParams, poly::DCRTPoly, sampler::uniform::DCRTPolyUniformSampler,
            FinRingElem,
        },
        sampler::DistType,
        PolyParams,
    };

    // Helper function to create a random polynomial using UniformSampler
    fn create_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
        let sampler = DCRTPolyUniformSampler::new();
        sampler.sample_poly(params, &DistType::FinRingDist)
    }

    fn create_bit_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
        let sampler = DCRTPolyUniformSampler::new();
        sampler.sample_poly(params, &DistType::BitDist)
    }

    #[test]
    fn test_eval_poly_circuit_add() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with an Add operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );

        // Expected result: poly1 + poly2
        let expected = poly1 + poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_eval_poly_circuit_sub() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with a Sub operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );

        // Expected result: poly1 - poly2
        let expected = poly1 - poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_eval_poly_circuit_scalar_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomial using UniformSampler
        let poly = create_random_poly(&params);

        // Create scalar as a random polynomial
        let scalar = create_random_poly(&params);

        // Create a circuit with a ScalarMul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(1);
        let scalar_mul_gate = circuit.scalar_mul_gate(inputs[0], scalar.clone());
        circuit.output(vec![scalar_mul_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &[poly.clone()]);

        // Expected result: poly * scalar
        let expected = poly * scalar;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_eval_poly_circuit_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with a Mul operation
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );

        // Expected result: poly1 * poly2
        let expected = poly1 * poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_eval_poly_circuit_complex() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);
        let poly3 = create_random_poly(&params);

        // Create scalar as a random polynomial
        let scalar = create_random_poly(&params);

        // Create a complex circuit: ((poly1 + poly2) * scalar) - poly3
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(3);

        // poly1 + poly2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (poly1 + poly2) * scalar
        let scalar_mul_gate = circuit.scalar_mul_gate(add_gate, scalar.clone());

        // ((poly1 + poly2) * scalar) - poly3
        let sub_gate = circuit.sub_gate(scalar_mul_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone(), poly3.clone()],
        );

        // Expected result: ((poly1 + poly2) * scalar) - poly3
        let expected = ((poly1 + poly2) * scalar) - poly3;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_eval_poly_circuit_multiple_outputs() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with multiple outputs
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);

        // poly1 + poly2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // poly1 - poly2
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);

        // poly1 * poly2
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);

        circuit.output(vec![add_gate, sub_gate, mul_gate]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );

        // Expected results
        let expected_add = poly1.clone() + poly2.clone();
        let expected_sub = poly1.clone() - poly2.clone();
        let expected_mul = poly1 * poly2;

        // Verify the results
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], expected_add);
        assert_eq!(result[1], expected_sub);
        assert_eq!(result[2], expected_mul);
    }

    #[test]
    fn test_eval_poly_circuit_deep_complex() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);
        let poly3 = create_random_poly(&params);
        let poly4 = create_random_poly(&params);

        // Create scalar as a random polynomial
        let scalar = create_random_poly(&params);

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = poly1 + poly2, b = poly3 * poly4
        // Level 2: c = a * b, d = poly1 - poly3
        // Level 3: e = c + d
        // Level 4: f = e * scalar
        // Output: f
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(4);

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // poly1 + poly2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // poly3 * poly4

        // Level 2
        let c = circuit.mul_gate(a, b); // (poly1 + poly2) * (poly3 * poly4)
        let d = circuit.sub_gate(inputs[0], inputs[2]); // poly1 - poly3

        // Level 3
        let e = circuit.add_gate(c, d); // ((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3)

        // Level 4
        let f = circuit.scalar_mul_gate(e, scalar.clone()); // (((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3)) * scalar

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone(), poly3.clone(), poly4.clone()],
        );

        // Expected result: (((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3)) * scalar
        let expected =
            (((poly1.clone() + poly2) * (poly3.clone() * poly4)) + (poly1 - poly3)) * scalar;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_boolean_gate_and() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let and_result = circuit.and_gate(inputs[0], inputs[1]);
        circuit.output(vec![and_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );
        let expected = poly1.clone() * poly2;
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_not() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(1);
        let not_result = circuit.not_gate(inputs[0]);
        circuit.output(vec![not_result]);
        let poly1 = create_bit_random_poly(&params);
        let result =
            circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &[poly1.clone()]);
        let expected = DCRTPoly::const_one(&params) - poly1.clone();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_or() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let or_result = circuit.or_gate(inputs[0], inputs[1]);
        circuit.output(vec![or_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );
        let expected = (poly1.clone() + poly2.clone()) - (poly1 * poly2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_nand() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let nand_result = circuit.nand_gate(inputs[0], inputs[1]);
        circuit.output(vec![nand_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );
        let expected = DCRTPoly::const_one(&params) - (poly1 * poly2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_nor() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let nor_result = circuit.nor_gate(inputs[0], inputs[1]); // poly1 AND poly2
        circuit.output(vec![nor_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );
        let expected =
            DCRTPoly::const_one(&params) - ((poly1.clone() + poly2.clone()) - (poly1 * poly2));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_xor() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let nor_result = circuit.xor_gate(inputs[0], inputs[1]);
        circuit.output(vec![nor_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );
        let expected = (poly1.clone() + poly2.clone())
            - (DCRTPoly::from_const(&params, &FinRingElem::new(2, params.modulus()))
                * poly1
                * poly2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_xnor() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(2);
        let xnor_result = circuit.xnor_gate(inputs[0], inputs[1]);
        circuit.output(vec![xnor_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result = circuit.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );
        let expected = DCRTPoly::const_one(&params)
            - ((poly1.clone() + poly2.clone())
                - (DCRTPoly::from_const(&params, &FinRingElem::new(2, params.modulus()))
                    * poly1
                    * poly2));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_mul_poly_bits_by_scalar_circuit() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let m = create_random_poly(&params);
        let decomposed_m = m.decompose(&params);
        let inputs = circuit.input(decomposed_m.len());
        assert_eq!(inputs.len(), params.modulus_bits());
        let k = create_bit_random_poly(&params); // Scalar

        // Output: decomposed_m[0] * k, ..., decomposed_m[modulus_bits - 1] * k
        let output_ids =
            inputs.iter().map(|&input_id| circuit.scalar_mul_gate(input_id, k.clone())).collect();

        circuit.output(output_ids);

        let result =
            circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &decomposed_m);

        assert_eq!(result.len(), params.modulus_bits());

        let recomposed_result = DCRTPoly::from_decomposed(&params, &result);
        assert_eq!(recomposed_result, m * k);
    }

    #[test]
    fn test_mul_fhe_poly_bits_mul_by_poly_circuit() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let uniform_sampler = DCRTPolyUniformSampler::new();

        // encrypt a polynomial m using a RLWE secret key encryption
        // c0 = a*s + e + m (where m is the plaintext polynomial)
        // c1 = -a
        let m = create_bit_random_poly(&params);
        let s = uniform_sampler.sample_poly(&params, &DistType::BitDist);
        let e = uniform_sampler.sample_poly(&params, &DistType::GaussDist { sigma: 0.0 });
        let a = uniform_sampler.sample_poly(&params, &DistType::FinRingDist);
        let c0 = -a.clone();
        let c1 = a * s.clone() + e + m.clone();

        // k is a polynomial from bit distribution
        let k = uniform_sampler.sample_poly(&params, &DistType::BitDist);

        let c0_bits = c0.decompose(&params);
        let c1_bits = c1.decompose(&params);

        let inputs = circuit.input(c0_bits.len() + c1_bits.len() + 1);
        assert_eq!(inputs.len(), params.modulus_bits() * 2 + 1);

        // Input: c0_bits[0], ..., c0_bits[modulus_bits - 1], c1_bits[0], ..., c1_bits[modulus_bits - 1], k
        // Output: c0_bits[0] * k, ..., c0_bits[modulus_bits - 1] * k, c1_bits[0] * k, ..., c1_bits[modulus_bits - 1] * k
        let k_id = inputs[inputs.len() - 1];
        let output_ids = inputs
            .iter()
            .take(inputs.len() - 1)
            .map(|&input_id| circuit.mul_gate(input_id, k_id))
            .collect();

        circuit.output(output_ids);

        // concatenate decomposed_c0 and decomposed_c1 and k
        let input = [c0_bits, c1_bits, vec![k.clone()]].concat();
        let result = circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &input);

        assert_eq!(result.len(), params.modulus_bits() * 2);

        let c0_bits_eval = result[..params.modulus_bits()].to_vec();
        let c1_bits_eval = result[params.modulus_bits()..].to_vec();

        let c0_eval = DCRTPoly::from_decomposed(&params, &c0_bits_eval);
        let c1_eval = DCRTPoly::from_decomposed(&params, &c1_bits_eval);

        assert_eq!(c0_eval, c0.clone() * k.clone());
        assert_eq!(c1_eval, c1.clone() * k.clone());

        // decrypt the result
        let plaintext = c1_eval + c0_eval * s;
        assert_eq!(plaintext, m * k);
    }
}
