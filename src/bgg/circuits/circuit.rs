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

        let mut input_gates = Vec::with_capacity(num_input);
        for idx in 0..num_input {
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

    pub fn add_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Add)
    }

    pub fn sub_gate(&mut self, left_input: usize, right_input: usize) -> usize {
        self.new_gate_generic(vec![left_input, right_input], PolyGateType::Sub)
    }

    pub fn scalar_mul_gate(&mut self, input: usize, scalar: P) -> usize {
        self.new_gate_generic(vec![input], PolyGateType::ScalarMul(scalar))
    }

    pub fn mul_gate(&mut self, left_input: usize, right_input: usize) -> usize {
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
}
