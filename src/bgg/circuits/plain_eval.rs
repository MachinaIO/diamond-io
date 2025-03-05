use itertools::Itertools;

use super::*;
use crate::poly::Poly;
use std::collections::BTreeMap;

pub fn eval_poly_circuit<P: Poly>(circuit: PolyCircuit<P>, input: &[P]) -> Vec<P> {
    #[cfg(debug_assertions)]
    {
        assert_eq!(circuit.num_input(), input.len());
    }

    let mut wires = BTreeMap::new();
    for (idx, input) in input.iter().enumerate() {
        wires.insert(idx, input.clone());
    }
    let output_gates = &circuit.output_ids.iter().map(|id| circuit.gates[id].clone()).collect_vec();
    for gate in output_gates.iter() {
        eval_poly_gate(&circuit, &mut wires, gate);
    }
    let output = circuit.output_ids.iter().map(|id| wires.get(id).unwrap().clone()).collect_vec();
    output
}

fn eval_poly_gate<P: Poly>(
    circuit: &PolyCircuit<P>,
    wires: &mut BTreeMap<usize, P>,
    gate: &PolyGate<P>,
) {
    let input_ids = &gate.input_gates;
    for input_id in input_ids.iter() {
        if !wires.contains_key(input_id) {
            let input_gate = &circuit.gates[input_id];
            eval_poly_gate(circuit, wires, input_gate);
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
            let output = input.clone() * scalar;
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
