use super::{PolyCircuit, PolyGateType};
use crate::poly::Poly;
use bytes::Bytes;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SerializablePolyGateType {
    Input,
    Add,
    Sub,
    ScalarMul(Bytes),
    Mul,
    Call { circuit_id: usize, num_input: usize, output_id: usize },
}

impl SerializablePolyGateType {
    pub fn num_input(&self) -> usize {
        match self {
            SerializablePolyGateType::Input => 0,
            SerializablePolyGateType::ScalarMul(_) => 1,
            SerializablePolyGateType::Add
            | SerializablePolyGateType::Sub
            | SerializablePolyGateType::Mul => 2,
            SerializablePolyGateType::Call { num_input, .. } => *num_input,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SerializablePolyGate {
    pub gate_id: usize,
    pub gate_type: SerializablePolyGateType,
    pub input_gates: Vec<usize>,
}

impl SerializablePolyGate {
    pub fn new(
        gate_id: usize,
        gate_type: SerializablePolyGateType,
        input_gates: Vec<usize>,
    ) -> Self {
        Self { gate_id, gate_type, input_gates }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SerializablePolyCircuit {
    gates: BTreeMap<usize, SerializablePolyGate>,
    sub_circuits: BTreeMap<usize, Self>,
    output_ids: Vec<usize>,
    num_input: usize,
}

impl SerializablePolyCircuit {
    pub fn new(
        gates: BTreeMap<usize, SerializablePolyGate>,
        sub_circuits: BTreeMap<usize, Self>,
        output_ids: Vec<usize>,
        num_input: usize,
    ) -> Self {
        Self { gates, sub_circuits, output_ids, num_input }
    }

    pub fn from_circuit<P: Poly>(circuit: &PolyCircuit<P>) -> Self {
        let mut gates = BTreeMap::new();
        for (gate_id, gate) in circuit.gates.iter() {
            let gate_type = match gate.gate_type {
                PolyGateType::Input => SerializablePolyGateType::Input,
                PolyGateType::Add => SerializablePolyGateType::Add,
                PolyGateType::Sub => SerializablePolyGateType::Sub,
                PolyGateType::ScalarMul(ref scalar) => SerializablePolyGateType::ScalarMul(
                    Bytes::from(scalar.to_compact_bytes(byte_size)),
                ),
            };
            let serializable_gate = SerializablePolyGate::new(
                *gate_id,
                SerializablePolyGateType::Input,
                gate.input_gates.clone(),
            );
            gates.insert(gate_id, serializable_gate);
        }

        let mut sub_circuits = BTreeMap::new();
        for (circuit_id, sub_circuit) in circuit.sub_circuits.iter() {
            let serializable_sub_circuit = Self::from_circuit(sub_circuit.clone());
            sub_circuits.insert(*circuit_id, serializable_sub_circuit);
        }
        Self::new(gates, sub_circuits, circuit.output_ids.clone(), circuit.num_input)
    }
}
