pub mod eval;
pub mod gate;
pub mod serde;
pub mod utils;
use dashmap::DashMap;
pub use eval::*;
pub use gate::{PolyGate, PolyGateType};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
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

    pub fn const_digits_poly(&mut self, digits: &[u32]) -> usize {
        self.new_gate_generic(vec![], PolyGateType::Const { digits: digits.to_vec() })
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

    /// Computes a topological order (as a vector of gate IDs) for all gates that
    /// are needed to evaluate the outputs. This is done via a DFS from each output gate.
    fn topological_order(&self) -> Vec<usize> {
        let mut visited = HashSet::new();
        let mut order = Vec::new();
        let mut stack = Vec::new();
        for &output_gate in &self.output_ids {
            if visited.insert(output_gate) {
                stack.push((output_gate, 0));
            }
        }

        while let Some((node, child_idx)) = stack.pop() {
            let gate = self.gates.get(&node).expect("gate not found");

            if child_idx < gate.input_gates.len() {
                stack.push((node, child_idx + 1));
                let child = gate.input_gates[child_idx];
                if visited.insert(child) {
                    stack.push((child, 0));
                }
            } else {
                order.push(node);
            }
        }

        order
    }

    /// Computes a levelized grouping of gate ids.
    /// Input wires (keys 0..=num_input) are assigned level 0.
    /// Each non‐input gate's level is defined as max(levels of its inputs) + 1.
    fn compute_levels(&self) -> Vec<Vec<usize>> {
        let mut gate_levels: HashMap<usize, usize> = HashMap::new();
        let mut levels: Vec<Vec<usize>> = vec![vec![]];
        for i in 0..=self.num_input {
            gate_levels.insert(i, 0);
        }
        let orders = self.topological_order();
        for gate_id in orders.into_iter() {
            let gate = self.gates.get(&gate_id).expect("gate not found");
            if gate.input_gates.is_empty() {
                gate_levels.insert(gate_id, 0);
                levels[0].push(gate_id);
                continue;
            }
            let level = gate
                .input_gates
                .iter()
                .map(|id| *gate_levels.get(id).expect("input gate not found"))
                .max()
                .expect("max level not found") +
                1;
            gate_levels.insert(gate_id, level);
            if levels.len() <= level {
                levels.resize(level + 1, vec![]);
            }
            levels[level].push(gate_id);
        }
        levels
    }

    /// Evaluate the circuit using an iterative approach over a precomputed topological order.
    pub fn eval<E: Evaluable>(&self, params: &E::Params, one: &E, inputs: &[E]) -> Vec<E> {
        #[cfg(debug_assertions)]
        {
            assert_eq!(self.num_input(), inputs.len());
            assert_ne!(self.num_output(), 0);
        }

        let wires = DashMap::new();
        let levels = self.compute_levels();
        debug_mem(format!("Levels: {:?}", levels));
        debug_mem("Levels are computed");

        wires.insert(0, one.clone());
        for (idx, input) in inputs.iter().enumerate() {
            wires.insert(idx + 1, input.clone());
        }
        debug_mem("Input wires are set");

        for level in levels.iter() {
            debug_mem("New level started");
            // All gates in the same level can be processed in parallel.
            level.par_iter().for_each(|&gate_id| {
                debug_mem(format!("Gate id {} started", gate_id));
                if wires.contains_key(&gate_id) {
                    debug_mem(format!("Gate id {} already evaluated", gate_id));
                    return;
                }
                let gate = self.gates.get(&gate_id).expect("gate not found").clone();
                debug_mem("Get gate");
                let result = match &gate.gate_type {
                    PolyGateType::Input => {
                        panic!("Input gate {:?} should already be preloaded", gate);
                    }
                    PolyGateType::Const { digits } => E::from_digits(params, one, digits),
                    PolyGateType::Add => {
                        debug_mem("Add gate start");
                        let left =
                            wires.get(&gate.input_gates[0]).expect("wire missing for Add").clone();
                        let right =
                            wires.get(&gate.input_gates[1]).expect("wire missing for Add").clone();
                        let result = left + right;
                        debug_mem("Add gate end");
                        result
                    }
                    PolyGateType::Sub => {
                        debug_mem("Sub gate start");
                        let left =
                            wires.get(&gate.input_gates[0]).expect("wire missing for Sub").clone();
                        let right =
                            wires.get(&gate.input_gates[1]).expect("wire missing for Sub").clone();
                        let result = left - right;
                        debug_mem("Sub gate end");
                        result
                    }
                    PolyGateType::Mul => {
                        debug_mem("Mul gate start");
                        let left =
                            wires.get(&gate.input_gates[0]).expect("wire missing for Mul").clone();
                        let right =
                            wires.get(&gate.input_gates[1]).expect("wire missing for Mul").clone();
                        let result = left * right;
                        debug_mem("Mul gate end");
                        result
                    }
                    PolyGateType::Rotate { shift } => {
                        debug_mem("Rotate gate start");
                        let input =
                            wires.get(&gate.input_gates[0]).expect("wire missing for Rotate");
                        let result = input.rotate(params, *shift);
                        debug_mem("Rotate gate end");
                        result
                    }
                    PolyGateType::Call { .. } => {
                        panic!("no more call gate type during evaluation");
                    }
                };
                wires.insert(gate_id, result);
                debug_mem(format!("Gate id {} finished", gate_id));
            });
            debug_mem("Evaluated gate in parallel");
        }

        let outputs = self
            .output_ids
            .par_iter()
            .map(|&id| wires.get(&id).expect("output missing").clone())
            .collect();
        debug_mem("Outputs are collected");
        outputs
    }

    pub fn register_sub_circuit(&mut self, sub_circuit: Self) -> usize {
        let circuit_id = self.sub_circuits.len();
        self.sub_circuits.insert(circuit_id, sub_circuit);
        circuit_id
    }

    /// Inlines the subcircuit operations directly into the main circuit instead of using call
    /// gates.
    pub fn call_sub_circuit(&mut self, circuit_id: usize, inputs: &[usize]) -> Vec<usize> {
        #[cfg(debug_assertions)]
        {
            let sub_circuit = &self.sub_circuits[&circuit_id];
            assert_eq!(inputs.len(), sub_circuit.num_input());
        }
        let mut gate_map: BTreeMap<usize, usize> = BTreeMap::new();
        let sub_circuit = self.sub_circuits.get(&circuit_id).unwrap().clone();
        for i in 0..=sub_circuit.num_input {
            if i == 0 {
                gate_map.insert(i, 0);
            } else if i <= inputs.len() {
                gate_map.insert(i, inputs[i - 1]);
            }
        }

        let mut outputs = Vec::with_capacity(sub_circuit.num_output());
        for &output_id in &sub_circuit.output_ids {
            let main_gate_id = self.inline_gate(output_id, &sub_circuit, &mut gate_map);
            outputs.push(main_gate_id);
        }
        outputs
    }

    /// Iteratively inlines a gate and its dependencies from a subcircuit into the main circuit.
    /// Returns the ID of the corresponding gate in the main circuit.
    fn inline_gate(
        &mut self,
        start_gate_id: usize,
        sub_circuit: &PolyCircuit,
        gate_map: &mut BTreeMap<usize, usize>,
    ) -> usize {
        let mut stack = Vec::new();
        stack.push(start_gate_id);

        while let Some(&current_gate_id) = stack.last() {
            if gate_map.contains_key(&current_gate_id) {
                stack.pop();
                continue;
            }
            let gate = sub_circuit.gates.get(&current_gate_id).unwrap();
            let mut all_inputs_inlined = true;
            for &input_id in &gate.input_gates {
                if !gate_map.contains_key(&input_id) {
                    all_inputs_inlined = false;
                    stack.push(input_id);
                }
            }
            if all_inputs_inlined {
                let main_inputs: Vec<usize> =
                    gate.input_gates.iter().map(|input_id| gate_map[input_id]).collect();
                let main_gate_id = self.new_gate_generic(main_inputs, gate.gate_type.clone());
                gate_map.insert(current_gate_id, main_gate_id);
                stack.pop();
            }
        }
        gate_map[&start_gate_id]
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
                DCRTPolyMatrix, FinRingElem,
            },
            enc::rlwe_encrypt,
            sampler::DistType,
            Poly, PolyMatrix, PolyParams,
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
    fn test_const_digits_poly() {
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
        let digits = vec![1u32, 0u32, 1u32, 1u32];
        let digits_poly_gate = circuit.const_digits_poly(&digits);
        circuit.output(vec![digits_poly_gate]);

        // Evaluate the circuit with any input (it won't be used)
        let dummy_input = create_random_poly(&params);
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &[dummy_input]);

        // Verify the result
        assert_eq!(result.len(), 1);

        // Check that the coefficients match the bit pattern
        let coeffs = result[0].coeffs();
        for (i, digit) in digits.iter().enumerate() {
            if digit != &0 {
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
        for (i, _) in coeffs.iter().enumerate().skip(digits.len()) {
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
        let mut circuit = PolyCircuit::new();
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();
        let sigma = 3.0;
        let log_q = params.modulus_bits();

        // encrypt a polynomial m using RLWE secret key encryption
        // Generate random message bits
        let m = sampler.sample_poly(&params, &DistType::BitDist);

        // Encrypt the message
        let a = sampler.sample_poly(&params, &DistType::BitDist);
        let t = sampler.sample_poly(&params, &DistType::BitDist);

        let m_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![m.clone()]);
        let a_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![a.clone()]);
        let t_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![t.clone()]);
        let b_mat = rlwe_encrypt(&params, &sampler, &t_mat, &a_mat, &m_mat, sigma);
        let b = b_mat.entry(0, 0);

        // ct = (a, b)
        let a_bits = a.decompose_base(&params);
        let b_bits = b.decompose_base(&params);

        let x = DCRTPoly::const_one(&params);

        let inputs = circuit.input(a_bits.len() + b_bits.len() + 1);
        assert_eq!(inputs.len(), params.modulus_bits() * 2 + 1);

        // Input: ct[bits], x
        // Output: ct[bits] * x
        let x_id = inputs[inputs.len() - 1];
        let output_ids = inputs
            .iter()
            .take(inputs.len() - 1)
            .map(|&input_id| circuit.mul_gate(input_id, x_id))
            .collect();

        circuit.output(output_ids);

        // concatenate decomposed_c0 and decomposed_c1 and x
        let input = [a_bits, b_bits, vec![x.clone()]].concat();
        let result = circuit.eval(&params, &DCRTPoly::const_one(&params), &input);

        assert_eq!(result.len(), log_q * 2);

        let a_bits_eval = result[..params.modulus_bits()].to_vec();
        let b_bits_eval = result[params.modulus_bits()..].to_vec();

        let a_eval = DCRTPoly::from_decomposed(&params, &a_bits_eval);
        let b_eval = DCRTPoly::from_decomposed(&params, &b_bits_eval);

        assert_eq!(a_eval, &a * &x);
        assert_eq!(b_eval, &b * &x);

        // decrypt the result
        let plaintext = b_eval - (a_eval * t);
        // recover the message bits
        let plaintext_bits = plaintext.extract_bits_with_threshold(&params);
        assert_eq!(plaintext_bits, (m * x).to_bool_vec());
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
