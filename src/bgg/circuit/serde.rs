use super::{PolyCircuit, PolyGateType};
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SerializablePolyGateType {
    Input,
    Const { bits: Vec<bool> },
    Add,
    Sub,
    Mul,
    Rotate { shift: usize },
    Call { circuit_id: usize, num_input: usize, output_id: usize },
}

impl SerializablePolyGateType {
    pub fn num_input(&self) -> usize {
        match self {
            SerializablePolyGateType::Input | SerializablePolyGateType::Const { .. } => 0,
            SerializablePolyGateType::Rotate { .. } => 1,
            SerializablePolyGateType::Add |
            SerializablePolyGateType::Sub |
            SerializablePolyGateType::Mul => 2,
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

    pub fn from_circuit(circuit: &PolyCircuit) -> Self {
        let mut gates = BTreeMap::new();
        for (gate_id, gate) in circuit.gates.iter() {
            let gate_type = match &gate.gate_type {
                PolyGateType::Input => SerializablePolyGateType::Input,
                PolyGateType::Const { bits } => {
                    SerializablePolyGateType::Const { bits: bits.clone() }
                }
                PolyGateType::Add => SerializablePolyGateType::Add,
                PolyGateType::Sub => SerializablePolyGateType::Sub,
                PolyGateType::Mul => SerializablePolyGateType::Mul,
                PolyGateType::Rotate { shift } => {
                    SerializablePolyGateType::Rotate { shift: *shift }
                }
                PolyGateType::Call { circuit_id, num_input, output_id } => {
                    SerializablePolyGateType::Call {
                        circuit_id: *circuit_id,
                        num_input: *num_input,
                        output_id: *output_id,
                    }
                }
            };
            let serializable_gate =
                SerializablePolyGate::new(*gate_id, gate_type, gate.input_gates.clone());
            gates.insert(*gate_id, serializable_gate);
        }

        let mut sub_circuits = BTreeMap::new();
        for (circuit_id, sub_circuit) in circuit.sub_circuits.iter() {
            let serializable_sub_circuit = Self::from_circuit(&sub_circuit);
            sub_circuits.insert(*circuit_id, serializable_sub_circuit);
        }
        Self::new(gates, sub_circuits, circuit.output_ids.clone(), circuit.num_input)
    }

    pub fn to_circuit(&self) -> PolyCircuit {
        let mut circuit = PolyCircuit::new();
        circuit.input(self.num_input);
        for (_, serializable_sub_circuit) in self.sub_circuits.iter() {
            let sub_circuit = serializable_sub_circuit.to_circuit();
            circuit.register_sub_circuit(sub_circuit);
        }

        // Process gates in ascending order of their usize keys
        let mut gate_idx = 0;
        while gate_idx < self.gates.len() {
            let serializable_gate = &self.gates[&gate_idx];
            match &serializable_gate.gate_type {
                SerializablePolyGateType::Input => {
                    gate_idx += 1;
                }
                SerializablePolyGateType::Const { bits } => {
                    circuit.const_bit_poly(&bits);
                    gate_idx += 1;
                }
                SerializablePolyGateType::Add => {
                    circuit.add_gate(
                        serializable_gate.input_gates[0],
                        serializable_gate.input_gates[1],
                    );
                    gate_idx += 1;
                }
                SerializablePolyGateType::Sub => {
                    circuit.sub_gate(
                        serializable_gate.input_gates[0],
                        serializable_gate.input_gates[1],
                    );
                    gate_idx += 1;
                }
                SerializablePolyGateType::Mul => {
                    circuit.mul_gate(
                        serializable_gate.input_gates[0],
                        serializable_gate.input_gates[1],
                    );
                    gate_idx += 1;
                }
                SerializablePolyGateType::Rotate { shift } => {
                    circuit.rotate_gate(serializable_gate.input_gates[0], *shift);
                    gate_idx += 1;
                }
                SerializablePolyGateType::Call { circuit_id, .. } => {
                    let output_size = circuit.sub_circuits[circuit_id].num_output();
                    circuit.call_sub_circuit(*circuit_id, &serializable_gate.input_gates);
                    gate_idx += output_size;
                }
            };
        }
        circuit.output(self.output_ids.clone());
        circuit
    }

    pub fn from_json_str(json_str: &str) -> Self {
        serde_json::from_str(json_str).expect("Failed to deserialize SerializablePolyCircuit")
    }

    pub fn to_json_str(&self) -> String {
        serde_json::to_string(self).expect("Failed to serialize SerializablePolyCircuit")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serialization_roundtrip() {
        // Create a complex circuit with various operations
        let mut original_circuit = PolyCircuit::new();

        // Add inputs
        let inputs = original_circuit.input(3);

        // Add various gates
        let add_gate = original_circuit.add_gate(inputs[0], inputs[1]);
        let sub_gate = original_circuit.sub_gate(add_gate, inputs[2]);
        let mul_gate = original_circuit.mul_gate(inputs[1], inputs[2]);

        // Create a sub-circuit
        let mut sub_circuit = PolyCircuit::new();
        let sub_inputs = sub_circuit.input(2);
        let sub_add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);
        let sub_mul_gate = sub_circuit.mul_gate(sub_inputs[0], sub_inputs[1]);
        sub_circuit.output(vec![sub_add_gate, sub_mul_gate]);

        // Register the sub-circuit
        let sub_circuit_id = original_circuit.register_sub_circuit(sub_circuit);

        // Call the sub-circuit
        let sub_outputs =
            original_circuit.call_sub_circuit(sub_circuit_id, &[inputs[0], inputs[1]]);

        // Use the sub-circuit outputs
        let combined_gate = original_circuit.add_gate(sub_gate, sub_outputs[0]);

        // Set the output
        original_circuit.output(vec![combined_gate, mul_gate, sub_outputs[1]]);

        // Convert to SerializablePolyCircuit
        let serializable_circuit = SerializablePolyCircuit::from_circuit(&original_circuit);

        // Convert back to PolyCircuit
        let roundtrip_circuit = serializable_circuit.to_circuit();

        // Verify that the circuits are identical by directly comparing them
        // This works because PolyCircuit implements the Eq trait
        assert_eq!(roundtrip_circuit, original_circuit);
    }

    #[test]
    fn test_serialization_roundtrip_json() {
        // Create a complex circuit with various operations
        let mut original_circuit = PolyCircuit::new();

        // Add inputs
        let inputs = original_circuit.input(3);

        // Add various gates
        let add_gate = original_circuit.add_gate(inputs[0], inputs[1]);
        let mul_gate = original_circuit.mul_gate(inputs[1], inputs[2]);

        // Create a sub-circuit
        let mut sub_circuit = PolyCircuit::new();
        let sub_inputs = sub_circuit.input(2);
        let sub_add_gate = sub_circuit.add_gate(sub_inputs[0], sub_inputs[1]);
        let sub_mul_gate = sub_circuit.mul_gate(sub_inputs[0], sub_inputs[1]);
        sub_circuit.output(vec![sub_add_gate, sub_mul_gate]);

        // Register the sub-circuit
        let sub_circuit_id = original_circuit.register_sub_circuit(sub_circuit);

        // Call the sub-circuit
        let sub_outputs =
            original_circuit.call_sub_circuit(sub_circuit_id, &[inputs[0], inputs[1]]);

        // Use the sub-circuit outputs
        let combined_gate = original_circuit.add_gate(add_gate, sub_outputs[0]);

        // Set the output
        original_circuit.output(vec![combined_gate, mul_gate, sub_outputs[1]]);

        // Convert to SerializablePolyCircuit
        let serializable_circuit = SerializablePolyCircuit::from_circuit(&original_circuit);
        let serializable_circuit_json = serializable_circuit.to_json_str();
        println!("{}", serializable_circuit_json);
        // Convert back to PolyCircuit
        let serializable_circuit =
            SerializablePolyCircuit::from_json_str(&serializable_circuit_json);
        let roundtrip_circuit = serializable_circuit.to_circuit();

        // Verify that the circuits are identical by directly comparing them
        // This works because PolyCircuit implements the Eq trait
        assert_eq!(roundtrip_circuit, original_circuit);
    }
}
