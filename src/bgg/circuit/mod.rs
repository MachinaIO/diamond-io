pub mod eval;
pub mod gate;
pub mod serde;
pub mod utils;
use dashmap::DashMap;
pub use eval::*;
pub use gate::{PolyGate, PolyGateType};
use itertools::Itertools;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    collections::{BTreeMap, HashMap},
    fmt::Debug,
};
pub use utils::*;

use crate::utils::debug_mem;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PolyCircuit {
    gates: BTreeMap<usize, PolyGate>,
    sub_circuits: BTreeMap<usize, Self>,
    output_ids: Vec<usize>,
    num_input: usize,
}

impl PolyCircuit {
    pub fn new() -> Self {
        Self {
            gates: BTreeMap::new(),
            sub_circuits: BTreeMap::new(),
            output_ids: vec![],
            num_input: 0,
        }
    }

    /// Get number of inputs
    pub fn num_input(&self) -> usize {
        self.num_input
    }

    /// Get number of outputs
    pub fn num_output(&self) -> usize {
        self.output_ids.len()
    }

    /// Get number of gates
    pub fn num_gates(&self) -> usize {
        self.gates.len()
    }

    pub fn input(&mut self, num_input: usize) -> Vec<usize> {
        #[cfg(debug_assertions)]
        assert_eq!(self.num_input, 0);
        self.gates.insert(0, PolyGate::new(0, PolyGateType::Input, vec![])); // input gate at index 0 reserved for constant 1 polynomial
        let mut input_gates = Vec::with_capacity(num_input);
        for idx in 1..(num_input + 1) {
            self.gates.insert(idx, PolyGate::new(idx, PolyGateType::Input, vec![]));
            input_gates.push(idx);
        }
        self.num_input = num_input;
        input_gates
    }

    pub fn output(&mut self, outputs: Vec<usize>) {
        #[cfg(debug_assertions)]
        assert_eq!(self.output_ids.len(), 0);

        for gate_id in outputs.into_iter() {
            self.output_ids.push(gate_id);
        }
    }

    pub fn const_zero_gate(&mut self) -> usize {
        self.not_gate(0)
    }

    /// index 0 have value 1
    pub fn const_one_gate(&mut self) -> usize {
        0
    }

    pub fn const_minus_one_gate(&mut self) -> usize {
        let zero = self.const_zero_gate();
        self.sub_gate(zero, 0)
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

    pub fn register_sub_circuit(&mut self, sub_circuit: Self) -> usize {
        let circuit_id = self.sub_circuits.len();
        self.sub_circuits.insert(circuit_id, sub_circuit);
        circuit_id
    }

    pub fn call_sub_circuit(&mut self, circuit_id: usize, inputs: &[usize]) -> Vec<usize> {
        let sub_circuit = &self.sub_circuits[&circuit_id];
        #[cfg(debug_assertions)]
        assert_eq!(inputs.len(), sub_circuit.num_input());
        let mut outputs = Vec::with_capacity(sub_circuit.num_output());
        let base_gate_id = self.gates.len();
        for idx in 0..sub_circuit.num_output() {
            let gate_id = self.new_gate_generic(
                inputs.to_vec(),
                PolyGateType::Call {
                    circuit_id,
                    num_input: inputs.len(),
                    output_id: base_gate_id + idx,
                },
            );
            outputs.push(gate_id);
        }
        outputs
    }

    pub fn add_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Add)
    }

    pub fn sub_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Sub)
    }

    pub fn mul_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Mul)
    }

    pub fn rotate_gate(&mut self, input: usize, shift: usize) -> usize {
        self.new_gate_generic(vec![input], PolyGateType::Rotate { shift })
    }

    pub fn const_bit_poly(&mut self, bits: &[bool]) -> usize {
        // let zero = self.const_zero_gate();
        // let one = self.const_one_gate();
        // let mut sum = zero.clone();
        // for (idx, bit_coeff) in coeffs.iter().enumerate() {
        //     let bit = if *bit_coeff { one } else { zero };
        //     let rotated = self.rotate_gate(bit, idx);
        //     sum = self.add_gate(sum, rotated);
        // }
        // sum
        self.new_gate_generic(vec![], PolyGateType::Const { bits: bits.to_vec() })
    }

    fn new_gate_generic(&mut self, inputs: Vec<usize>, gate_type: PolyGateType) -> usize {
        #[cfg(debug_assertions)]
        {
            assert_ne!(self.num_input, 0);
            assert_eq!(self.output_ids.len(), 0);
            assert_eq!(inputs.len(), gate_type.num_input());
            for gate_id in inputs.iter() {
                assert!(self.gates.contains_key(gate_id));
            }
        }
        let gate_id = self.gates.len();
        self.gates.insert(gate_id, PolyGate::new(gate_id, gate_type, inputs));
        gate_id
    }

    /// Computes a levelized grouping of gate ids.
    /// Input wires (keys 0..=num_input) are assigned level 0.
    /// Each non‐input gate’s level is defined as max(levels of its inputs) + 1.
    pub fn compute_levels(&self) -> Vec<Vec<usize>> {
        let mut gate_levels: HashMap<usize, usize> = HashMap::new();
        let mut levels: Vec<Vec<usize>> = vec![];

        // Input wires (including constant one at key 0) are at level 0.
        for i in 0..=self.num_input {
            gate_levels.insert(i, 0);
        }

        // For each gate, compute its level
        for (&gate_id, gate) in self.gates.iter() {
            // skipping input gates
            if gate_id <= self.num_input {
                continue;
            }

            let level = gate
                .input_gates
                .iter()
                .map(|id| {
                    let r = *gate_levels.get(id).unwrap_or(&0);
                    println!("gate: {} / {:?}", r, gate_levels.get(id));
                    r
                })
                .max()
                .unwrap_or(0) +
                1;
            gate_levels.insert(gate_id, level);
            if levels.len() <= level {
                levels.resize(level + 1, vec![]);
            }
            println!("level :{}, gate_id: {}", level, gate_id);
            levels[level].push(gate_id);
        }
        levels
    }

    /// Each gate is evaluated assuming that all its input wires have already been computed.
    fn eval_gate<E: Evaluable>(
        &self,
        params: &E::Params,
        one: &E,
        wires: &DashMap<usize, E>,
        gate: &PolyGate,
    ) -> E {
        let input_ids = &gate.input_gates;
        match &gate.gate_type {
            PolyGateType::Input => {
                panic!("Input gate should have been initialized")
            }
            PolyGateType::Const { bits } => E::from_bits(params, one, bits),
            PolyGateType::Add => {
                let left = wires.get(&input_ids[0]).expect("missing wire").clone();
                let right = wires.get(&input_ids[1]).expect("missing wire").clone();
                left + right
            }
            PolyGateType::Sub => {
                let left = wires.get(&input_ids[0]).expect("missing wire").clone();
                let right = wires.get(&input_ids[1]).expect("missing wire").clone();
                left - right
            }
            PolyGateType::Mul => {
                let left = wires.get(&input_ids[0]).expect("missing wire").clone();
                let right = wires.get(&input_ids[1]).expect("missing wire").clone();
                left * right
            }
            PolyGateType::Rotate { shift } => {
                let input = wires.get(&input_ids[0]).expect("missing wire").clone();
                input.rotate(params, *shift)
            }
            PolyGateType::Call { circuit_id, output_id, .. } => {
                let sub_circuit = self.sub_circuits.get(circuit_id).expect("missing sub-circuit");
                let sub_inputs: Vec<E> = input_ids
                    .iter()
                    .map(|id| wires.get(id).expect("missing sub-input").clone())
                    .collect();
                let outputs = sub_circuit.eval(params, one, &sub_inputs);
                for (idx, output_wire) in outputs.into_iter().enumerate() {
                    wires.insert(output_id + idx, output_wire);
                }
                wires.get(output_id).expect("sub-circuit output missing").clone()
            }
        }
    }

    /// Evaluate the circuit in parallel using a levelized approach.
    pub fn eval<E: Evaluable>(&self, params: &E::Params, one: &E, inputs: &[E]) -> Vec<E> {
        #[cfg(debug_assertions)]
        {
            assert_eq!(self.num_input(), inputs.len());
            assert_ne!(self.num_output(), 0);
        }
        let wires = DashMap::new();
        wires.insert(0, one.clone());
        for (idx, input) in inputs.iter().enumerate() {
            wires.insert(idx + 1, input.clone());
        }
        let levels = self.compute_levels();
        println!("levels :{:?}", levels);
        // Evaluate gates level-by-level.
        for level in levels.iter() {
            // All gates in the same level can be processed in parallel.
            level.par_iter().for_each(|&gate_id| {
                let gate = self.gates.get(&gate_id).expect("gate not found").clone();
                let res = self.eval_gate(params, one, &wires, &gate);
                wires.insert(gate.gate_id, res);
                debug_mem("Evaluated gate in parallel");
            });
        }

        self.output_ids
            .iter()
            .map(|id| wires.get(id).expect("output missing").clone())
            .collect_vec()
    }
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::{
        poly::{
            dcrt::{
                params::DCRTPolyParams, poly::DCRTPoly, sampler::uniform::DCRTPolyUniformSampler,
                FinRingElem,
            },
            sampler::DistType,
            Poly, PolyParams,
        },
        utils::{create_bit_random_poly, create_random_poly},
    };
    use num_bigint::BigUint;

    #[test]
    fn test_eval_add() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with an Add operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        circuit.output(vec![add_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);

        // Expected result: poly1 + poly2
        let expected = poly1 + poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_eval_sub() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with a Sub operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);
        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);

        // Expected result: poly1 - poly2
        let expected = poly1 - poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_eval_mul() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with a Mul operation
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);
        circuit.output(vec![mul_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);

        // Expected result: poly1 * poly2
        let expected = poly1 * poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_const_bit_poly() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a circuit with a const_bit_poly gate
        let mut circuit = PolyCircuit::new();
        // We need to call input() to initialize the circuit
        circuit.input(1);

        // Define a specific bit pattern
        // This will create a polynomial with coefficients:
        // [1, 0, 1, 1]
        // (where 1 is at positions 0, 2, 3, and 4)
        let bits = vec![true, false, true, true];
        let bit_poly_gate = circuit.const_bit_poly(&bits);
        circuit.output(vec![bit_poly_gate]);

        // Evaluate the circuit with any input (it won't be used)
        let dummy_input = create_random_poly(&params);
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &[dummy_input]);

        // Verify the result
        assert_eq!(result.len(), 1);

        // Check that the coefficients match the bit pattern
        let coeffs = result[0].coeffs();
        for (i, bit) in bits.iter().enumerate() {
            if *bit {
                assert_eq!(
                    coeffs[i].value(),
                    &BigUint::from(1u8),
                    "Coefficient at position {} should be 1",
                    i
                );
            } else {
                assert_eq!(
                    coeffs[i].value(),
                    &BigUint::from(0u8),
                    "Coefficient at position {} should be 0",
                    i
                );
            }
        }

        // Check that remaining coefficients are 0
        for (i, _) in coeffs.iter().enumerate().skip(bits.len()) {
            assert_eq!(
                coeffs[i].value(),
                &BigUint::from(0u8),
                "Coefficient at position {} should be 0",
                i
            );
        }
    }

    #[test]
    fn test_eval_complex() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);
        let poly3 = create_random_poly(&params);

        // Create a complex circuit: (poly1 + poly2) - poly3
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(3);

        // poly1 + poly2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // (poly1 + poly2) - poly3
        let sub_gate = circuit.sub_gate(add_gate, inputs[2]);

        circuit.output(vec![sub_gate]);

        // Evaluate the circuit
        let result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone(), poly3.clone()],
        );

        // Expected result: (poly1 + poly2) - poly3
        let expected = (poly1 + poly2) - poly3;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_eval_multiple_outputs() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a circuit with multiple outputs
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);

        // poly1 + poly2
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);

        // poly1 - poly2
        let sub_gate = circuit.sub_gate(inputs[0], inputs[1]);

        // poly1 * poly2
        let mul_gate = circuit.mul_gate(inputs[0], inputs[1]);

        circuit.output(vec![add_gate, sub_gate, mul_gate]);

        // Evaluate the circuit
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);

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
    fn test_eval_deep_complex() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);
        let poly3 = create_random_poly(&params);
        let poly4 = create_random_poly(&params);

        // Create a complex circuit with depth = 4
        // Circuit structure:
        // Level 1: a = poly1 + poly2, b = poly3 * poly4, d = poly1 - poly3
        // Level 2: c = a * b
        // Level 3: e = c + d
        // Level 4: f = e * e
        // Output: f
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(4);

        // Level 1
        let a = circuit.add_gate(inputs[0], inputs[1]); // poly1 + poly2
        let b = circuit.mul_gate(inputs[2], inputs[3]); // poly3 * poly4
        let d = circuit.sub_gate(inputs[0], inputs[2]); // poly1 - poly3

        // Level 2
        let c = circuit.mul_gate(a, b); // (poly1 + poly2) * (poly3 * poly4)

        // Level 3
        let e = circuit.add_gate(c, d); // ((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3)

        // Level 4
        let f = circuit.mul_gate(e, e); // (((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3))^2

        circuit.output(vec![f]);

        // Evaluate the circuit
        let result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone(), poly3.clone(), poly4.clone()],
        );

        // Expected result: (((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3))^2
        let expected = (((poly1.clone() + poly2.clone()) * (poly3.clone() * poly4.clone())) +
            (poly1.clone() - poly3.clone())) *
            (((poly1.clone() + poly2) * (poly3.clone() * poly4)) + (poly1 - poly3));

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_boolean_gate_and() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let and_result = circuit.and_gate(inputs[0], inputs[1]);
        circuit.output(vec![and_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);
        let expected = poly1.clone() * poly2;
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_not() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(1);
        let not_result = circuit.not_gate(inputs[0]);
        circuit.output(vec![not_result]);
        let poly1 = create_bit_random_poly(&params);
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone()]);
        let expected = DCRTPoly::const_one(&params) - poly1.clone();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_or() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let or_result = circuit.or_gate(inputs[0], inputs[1]);
        circuit.output(vec![or_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);
        let expected = (poly1.clone() + poly2.clone()) - (poly1 * poly2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_nand() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let nand_result = circuit.nand_gate(inputs[0], inputs[1]);
        circuit.output(vec![nand_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);
        let expected = DCRTPoly::const_one(&params) - (poly1 * poly2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_nor() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let nor_result = circuit.nor_gate(inputs[0], inputs[1]); // poly1 AND poly2
        circuit.output(vec![nor_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);
        let expected =
            DCRTPoly::const_one(&params) - ((poly1.clone() + poly2.clone()) - (poly1 * poly2));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_xor() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let nor_result = circuit.xor_gate(inputs[0], inputs[1]);
        circuit.output(vec![nor_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);
        let expected = (poly1.clone() + poly2.clone()) -
            (DCRTPoly::from_const(&params, &FinRingElem::new(2, params.modulus())) *
                poly1 *
                poly2);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_boolean_gate_xnor() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
        let inputs = circuit.input(2);
        let xnor_result = circuit.xnor_gate(inputs[0], inputs[1]);
        circuit.output(vec![xnor_result]);
        let poly1 = create_bit_random_poly(&params);
        let poly2 = create_bit_random_poly(&params);
        let result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &[poly1.clone(), poly2.clone()]);
        let expected = DCRTPoly::const_one(&params) -
            ((poly1.clone() + poly2.clone()) -
                (DCRTPoly::from_const(&params, &FinRingElem::new(2, params.modulus())) *
                    poly1 *
                    poly2));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].coeffs(), expected.coeffs());
    }

    #[test]
    fn test_mul_fhe_poly_bits_mul_by_poly_circuit() {
        let params = DCRTPolyParams::default();
        let mut circuit = PolyCircuit::new();
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

        // Input: c0_bits[0], ..., c0_bits[modulus_bits - 1], c1_bits[0], ..., c1_bits[modulus_bits
        // - 1], k Output: c0_bits[0] * k, ..., c0_bits[modulus_bits - 1] * k, c1_bits[0] *
        // k, ..., c1_bits[modulus_bits - 1] * k
        let k_id = inputs[inputs.len() - 1];
        let output_ids = inputs
            .iter()
            .take(inputs.len() - 1)
            .map(|&input_id| circuit.mul_gate(input_id, k_id))
            .collect();

        circuit.output(output_ids);

        // concatenate decomposed_c0 and decomposed_c1 and k
        let input = [c0_bits, c1_bits, vec![k.clone()]].concat();
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &input);

        assert_eq!(result.len(), params.modulus_bits() * 2);

        let c0_bits_eval = result[..params.modulus_bits()].to_vec();
        let c1_bits_eval = result[params.modulus_bits()..].to_vec();

        let c0_eval = DCRTPoly::from_decomposed(&params, &c0_bits_eval);
        let c1_eval = DCRTPoly::from_decomposed(&params, &c1_bits_eval);

        assert_eq!(c0_eval, &c0 * &k);
        assert_eq!(c1_eval, &c1 * &k);

        // decrypt the result
        let plaintext = c1_eval + c0_eval * s;
        assert_eq!(plaintext, m * k);
    }

    #[test]
    fn test_register_and_call_sub_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials using UniformSampler
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);

        // Create a sub-circuit that performs addition and multiplication
        let mut sub_circuit = PolyCircuit::new();
        let sub_inputs = sub_circuit.input(2);

        // Add operation: poly1 + poly2
        let add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);

        // Mul operation: poly1 * poly2
        let mul_gate = sub_circuit.mul_gate(sub_inputs[0], sub_inputs[1]);

        // Set the outputs of the sub-circuit
        sub_circuit.output(vec![add_gate, mul_gate]);

        // Create the main circuit
        let mut main_circuit = PolyCircuit::new();
        let main_inputs = main_circuit.input(2);

        // Register the sub-circuit and get its ID
        let sub_circuit_id = main_circuit.register_sub_circuit(sub_circuit);

        // Call the sub-circuit with the main circuit's inputs
        let sub_outputs =
            main_circuit.call_sub_circuit(sub_circuit_id, &[main_inputs[0], main_inputs[1]]);

        // Verify we got two outputs from the sub-circuit
        assert_eq!(sub_outputs.len(), 2);

        // Use the sub-circuit outputs for further operations
        // For example, subtract the multiplication result from the addition result
        let final_gate = main_circuit.sub_gate(sub_outputs[0], sub_outputs[1]);

        // Set the output of the main circuit
        main_circuit.output(vec![final_gate]);

        // Evaluate the main circuit
        let result = main_circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone()],
        );

        // Expected result: (poly1 + poly2) - (poly1 * poly2)
        let expected = (poly1.clone() + poly2.clone()) - (poly1 * poly2);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_nested_sub_circuits() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials
        let poly1 = create_random_poly(&params);
        let poly2 = create_random_poly(&params);
        let poly3 = create_random_poly(&params);

        // Create the innermost sub-circuit that performs multiplication
        let mut inner_circuit = PolyCircuit::new();
        let inner_inputs = inner_circuit.input(2);
        let mul_gate = inner_circuit.mul_gate(inner_inputs[0], inner_inputs[1]);
        inner_circuit.output(vec![mul_gate]);

        // Create a middle sub-circuit that uses the inner sub-circuit
        let mut middle_circuit = PolyCircuit::new();
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
        let mut main_circuit = PolyCircuit::new();
        let main_inputs = main_circuit.input(3);

        // Register the middle circuit
        let middle_circuit_id = main_circuit.register_sub_circuit(middle_circuit);

        // Call the middle circuit with all inputs
        let middle_outputs = main_circuit
            .call_sub_circuit(middle_circuit_id, &[main_inputs[0], main_inputs[1], main_inputs[2]]);

        let scalar_mul_gate = main_circuit.mul_gate(middle_outputs[0], middle_outputs[0]);

        // Set the output of the main circuit
        main_circuit.output(vec![scalar_mul_gate]);

        // Evaluate the main circuit
        let result = main_circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[poly1.clone(), poly2.clone(), poly3.clone()],
        );

        // Expected result: ((poly1 * poly2) + poly3)^2
        let expected =
            ((poly1.clone() * poly2.clone()) + poly3.clone()) * ((poly1 * poly2) + poly3);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_const_zero_gate() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a circuit with a const_zero_gate
        let mut circuit = PolyCircuit::new();
        // We need to call input() to initialize the circuit
        circuit.input(1);
        let zero_gate = circuit.const_zero_gate();
        circuit.output(vec![zero_gate]);

        // Evaluate the circuit with any input (it won't be used)
        let dummy_input = create_random_poly(&params);
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &[dummy_input]);

        // Expected result: 0
        let expected = DCRTPoly::const_zero(&params);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_const_one_gate() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a circuit with a const_one_gate
        let mut circuit = PolyCircuit::new();
        // We need to call input() to initialize the circuit
        circuit.input(1);
        let one_gate = circuit.const_one_gate();
        circuit.output(vec![one_gate]);

        // Evaluate the circuit with any input (it won't be used)
        let dummy_input = create_random_poly(&params);
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &[dummy_input]);

        // Expected result: 1
        let expected = DCRTPoly::const_one(&params);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }

    #[test]
    fn test_const_minus_one_gate() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create a circuit with a const_minus_one_gate
        let mut circuit = PolyCircuit::new();
        // We need to call input() to initialize the circuit
        circuit.input(1);
        let minus_one_gate = circuit.const_minus_one_gate();
        circuit.output(vec![minus_one_gate]);

        // Evaluate the circuit with any input (it won't be used)
        let dummy_input = create_random_poly(&params);
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &[dummy_input]);

        // Expected result: -1
        // We can compute -1 as 0 - 1
        let expected = DCRTPoly::const_zero(&params) - DCRTPoly::const_one(&params);

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }
}
