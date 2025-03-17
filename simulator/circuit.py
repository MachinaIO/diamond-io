import json
from typing import Dict, List, Optional
import math
import unittest
import os


class CallData:
    def __init__(self, circuit_id: int, num_input: int, output_id: int):
        self.circuit_id = circuit_id
        self.num_input = num_input
        self.output_id = output_id


class SerializablePolyGateType:
    """Python equivalent of the Rust SerializablePolyGateType enum."""

    INPUT = "Input"
    ADD = "Add"
    SUB = "Sub"
    MUL = "Mul"
    ROTATE = "Rotate"
    CALL = "Call"


class NormBoundsPerGate:
    def __init__(self, h_norm: int, plaintext_norm: int):
        self.h_norm = h_norm
        self.plaintext_norm = plaintext_norm


class SerializablePolyGate:
    """Python equivalent of the Rust SerializablePolyGate struct."""

    def __init__(
        self,
        gate_id: int,
        gate_type: str,
        input_gates: List[int],
        rotate_data: Optional[int] = None,
        call_data: Optional[CallData] = None,
    ):
        """
        Initialize a SerializablePolyGate.

        Args:
            gate_id: The ID of the gate
            gate_type: The type of the gate (one of the SerializablePolyGateType constants)
            input_gates: List of input gate IDs
            rotate_data: Shift amount for Rotate gates
            call_data: Dictionary with circuit_id, num_input, and output_id for Call gates
        """
        self.gate_id = gate_id
        self.gate_type = gate_type
        self.input_gates = input_gates
        self.rotate_data = rotate_data
        self.call_data = call_data

    @classmethod
    def from_json(cls, json_data: Dict) -> "SerializablePolyGate":
        """
        Create a SerializablePolyGate from JSON data.

        Args:
            json_data: The JSON data to parse
        """
        gate_id = json_data["gate_id"]
        gate_type_data = json_data["gate_type"]
        input_gates = json_data["input_gates"]

        # Handle different gate types
        if isinstance(gate_type_data, dict):
            # Complex gate types are represented as objects
            if "Rotate" in gate_type_data:
                shift = gate_type_data["Rotate"]["shift"]

                return cls(
                    gate_id,
                    SerializablePolyGateType.ROTATE,
                    input_gates,
                    rotate_data=shift,
                )
            elif "Call" in gate_type_data:
                # Call gates have circuit_id, num_input, and output_id fields
                # call_data = gate_type_data["Call"]
                call_data = CallData(
                    gate_type_data["Call"]["circuit_id"],
                    gate_type_data["Call"]["num_input"],
                    gate_type_data["Call"]["output_id"],
                )
                return cls(
                    gate_id,
                    SerializablePolyGateType.CALL,
                    input_gates,
                    call_data=call_data,
                )
        else:
            # Simple gate types are represented as strings
            return cls(gate_id, gate_type_data, input_gates)


class SerializablePolyCircuit:
    """Python equivalent of the Rust SerializablePolyCircuit struct."""

    def __init__(
        self,
        gates: Dict[int, SerializablePolyGate],
        sub_circuits: Dict[int, "SerializablePolyCircuit"],
        output_ids: List[int],
        num_input: int,
        ring_dimension: int,
        modulus: int,
        d: int,
    ):
        """
        Initialize a SerializablePolyCircuit.

        Args:
            gates: Dictionary mapping gate IDs to SerializablePolyGate objects
            sub_circuits: Dictionary mapping circuit IDs to SerializablePolyCircuit objects
            output_ids: List of output gate IDs
            num_input: Number of input gates
            ring_dimension: The ring dimension for the polynomial
            modulus: The modulus for the polynomial
            d: The number of secret polynomials
        """
        self.gates = gates
        self.sub_circuits = sub_circuits
        self.output_ids = output_ids
        self.num_input = num_input
        self.ring_dimension = ring_dimension
        self.modulus = modulus
        self.d = d

    @classmethod
    def from_json(
        cls,
        json_data: Dict,
        ring_dimension: int,
        modulus: int,
        d: int,
    ) -> "SerializablePolyCircuit":
        """
        Create a SerializablePolyCircuit from JSON data.

        Args:
            json_data: The JSON data to parse
            ring_dimension: The ring dimension for the polynomial
            modulus: The modulus for the polynomial
        """
        # Parse gates
        gates = {}
        for gate_id_str, gate_data in json_data.get("gates", {}).items():
            gate_id = int(gate_id_str)
            gates[gate_id] = SerializablePolyGate.from_json(gate_data)

        # Parse sub-circuits recursively
        sub_circuits = {}
        for circuit_id_str, circuit_data in json_data.get("sub_circuits", {}).items():
            circuit_id = int(circuit_id_str)
            sub_circuits[circuit_id] = cls.from_json(
                circuit_data, ring_dimension, modulus, d
            )

        # Parse output IDs and num_input
        output_ids = json_data["output_ids"]
        num_input = json_data["num_input"]

        return cls(
            gates, sub_circuits, output_ids, num_input, ring_dimension, modulus, d
        )

    @classmethod
    def load_from_file(
        cls, file_path: str, ring_dimension: int, modulus: int, d: int
    ) -> "SerializablePolyCircuit":
        """
        Load a SerializablePolyCircuit from a JSON file.

        Args:
            file_path: Path to the JSON file
            ring_dimension: The ring dimension for the polynomial
            modulus: The modulus for the polynomial
        """
        with open(file_path, "r") as f:
            json_data = json.load(f)
        return cls.from_json(json_data, ring_dimension, modulus, d)

    def estimate_error_bound(self) -> List[NormBoundsPerGate]:
        """
        Estimate the error bound for the circuit.

        Returns:
            The error bound for the circuit
        """
        n = self.ring_dimension
        return self._estimate_error_bound_inner(
            NormBoundsPerGate(1, 1),
            [NormBoundsPerGate(1, n) for _ in range(self.num_input)],
        )

    def _estimate_error_bound_inner(
        self,
        one: NormBoundsPerGate,
        inputs: List[NormBoundsPerGate],
    ) -> List[NormBoundsPerGate]:
        n = self.ring_dimension
        log_q = math.ceil(math.log2(self.modulus))
        m = (self.d + 1) * log_q
        norms = {}
        norms[0] = one
        for i, input in enumerate(inputs):
            norms[i + 1] = input
        for gate_id in sorted(self.gates.keys()):
            gate = self.gates[gate_id]
            if gate.gate_type == SerializablePolyGateType.INPUT:
                continue
            elif gate.gate_type == SerializablePolyGateType.ADD:
                left = norms[gate.input_gates[0]]
                right = norms[gate.input_gates[1]]
                h_norm = left.h_norm + right.h_norm
                plaintext_norm = left.plaintext_norm + right.plaintext_norm
                norms[gate_id] = NormBoundsPerGate(h_norm, plaintext_norm)
            elif gate.gate_type == SerializablePolyGateType.SUB:
                left = norms[gate.input_gates[0]]
                right = norms[gate.input_gates[1]]
                h_norm = left.h_norm + right.h_norm
                plaintext_norm = left.plaintext_norm + right.plaintext_norm
                norms[gate_id] = NormBoundsPerGate(h_norm, plaintext_norm)
            elif gate.gate_type == SerializablePolyGateType.MUL:
                left = norms[gate.input_gates[0]]
                right = norms[gate.input_gates[1]]
                h_norm = left.h_norm * m * n + right.h_norm * left.plaintext_norm
                plaintext_norm = left.plaintext_norm * right.plaintext_norm
                norms[gate_id] = NormBoundsPerGate(h_norm, plaintext_norm)
            elif gate.gate_type == SerializablePolyGateType.ROTATE:
                input = norms[gate.input_gates[0]]
                norms[gate_id] = NormBoundsPerGate(input.h_norm, input.plaintext_norm)
            elif gate.gate_type == SerializablePolyGateType.CALL:
                circuit = self.sub_circuits[gate.call_data.circuit_id]
                sub_inputs = [norms[input_gate] for input_gate in gate.input_gates]
                outputs = circuit._estimate_error_bound_inner(one, sub_inputs)
                init_output_gate_id = gate_id - gate.call_data.output_id
                for i, output in enumerate(outputs):
                    output_gate_id = init_output_gate_id + i
                    norms[output_gate_id] = output
        return [norms[output_id] for output_id in self.output_ids]


class TestCircuit(unittest.TestCase):
    """Test cases for the SerializablePolyCircuit class."""

    def test_estimate_error_bound_with_circuit_test_json(self):
        """Test that estimate_error_bound works with circuit_test.json."""
        # Define parameters for the test
        ring_dimension = 1024
        modulus = 2**32
        d = 3  # Number of secret polynomials

        # Get the path to circuit_test.json
        current_dir = os.path.dirname(os.path.abspath(__file__))
        circuit_test_path = os.path.join(current_dir, "circuit_test.json")

        # Load the circuit from the JSON file
        circuit = SerializablePolyCircuit.load_from_file(
            circuit_test_path, ring_dimension, modulus, d
        )

        # Call estimate_error_bound and verify it runs without errors
        error_bounds = circuit.estimate_error_bound()
        print(error_bounds[0].h_norm)
        # Verify that the error bounds are valid
        self.assertIsNotNone(error_bounds)
        self.assertGreater(len(error_bounds), 0)

        # Check that each error bound has valid h_norm and plaintext_norm values
        for bound in error_bounds:
            self.assertIsInstance(bound, NormBoundsPerGate)
            self.assertGreaterEqual(bound.h_norm, 0)
            self.assertGreaterEqual(bound.plaintext_norm, 0)


if __name__ == "__main__":
    unittest.main()
