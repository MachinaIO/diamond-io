import json
from typing import Dict, List, Optional
from .from_compact_bytes import from_compact_bytes
import math


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
    SCALAR_MUL = "ScalarMul"
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
        scalar_mul_data: Optional[List[int]] = None,
        call_data: Optional[CallData] = None,
    ):
        """
        Initialize a SerializablePolyGate.

        Args:
            gate_id: The ID of the gate
            gate_type: The type of the gate (one of the SerializablePolyGateType constants)
            input_gates: List of input gate IDs
            scalar_mul_data: Coefficient array for ScalarMul gates (list of integers)
            call_data: Dictionary with circuit_id, num_input, and output_id for Call gates
        """
        self.gate_id = gate_id
        self.gate_type = gate_type
        self.input_gates = input_gates
        self.scalar_mul_data = scalar_mul_data
        self.call_data = call_data

    @classmethod
    def from_json(
        cls, json_data: Dict, ring_dimension: int, modulus: int
    ) -> "SerializablePolyGate":
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
            if "ScalarMul" in gate_type_data:
                # ScalarMul gates have an array of integers that need to be converted to bytes
                # and then processed with from_compact_bytes
                scalar_mul_array = gate_type_data["ScalarMul"]

                # Convert the array to bytes
                scalar_mul_bytes = bytes(scalar_mul_array)

                # Use from_compact_bytes to convert to coefficient array
                scalar_mul_data = from_compact_bytes(
                    scalar_mul_bytes, ring_dimension, modulus
                )

                return cls(
                    gate_id,
                    SerializablePolyGateType.SCALAR_MUL,
                    input_gates,
                    scalar_mul_data=scalar_mul_data,
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
            gates[gate_id] = SerializablePolyGate.from_json(
                gate_data, ring_dimension, modulus
            )

        # Parse sub-circuits recursively
        sub_circuits = {}
        for circuit_id_str, circuit_data in json_data.get("sub_circuits", {}).items():
            circuit_id = int(circuit_id_str)
            sub_circuits[circuit_id] = cls.from_json(
                circuit_data, ring_dimension, modulus
            )

        # Parse output IDs and num_input
        output_ids = json_data.get("output_ids", [])
        num_input = json_data.get("num_input", 0)

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

    @classmethod
    def estimate_error_bound(self) -> List[int]:
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

    @classmethod
    def _estimate_error_bound_inner(
        self,
        one: NormBoundsPerGate,
        inputs: List[NormBoundsPerGate],
        gadget_vector: List[int],
    ) -> List[int]:
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
            elif gate.gate_type == SerializablePolyGateType.SCALAR_MUL:
                # poly = norms[gate.input_gates[0]]
                # scalar = gate.scalar_mul_data[0]
                # h_norm = poly.h_norm * scalar
                # plaintext_norm = poly.plaintext_norm * scalar
                # norms[gate_id] = NormBoundsPerGate(h_norm, plaintext_norm)
                # [todo]
                pass
            elif gate.gate_type == SerializablePolyGateType.CALL:
                circuit: SerializablePolyCircuit = self.sub_circuits[
                    gate.call_data.circuit_id
                ]
                sub_inputs = [norms[input_gate] for input_gate in gate.input_gates]
                outputs = circuit._estimate_error_bound_inner(one, sub_inputs)
                init_output_gate_id = gate_id - gate.call_data.output_id
                for i, output in enumerate(outputs):
                    output_gate_id = init_output_gate_id + i
                    norms[output_gate_id] = output
        return [norms[output_id] for output_id in self.output_ids]

    # def evaluate(
    #     self, inputs: List[List[int]], modulus: int, ring_dimension: int
    # ) -> List[List[int]]:
    #     """
    #     Evaluate the circuit with the given inputs.

    #     This is a simplified evaluation function that only works for basic circuits.
    #     For more complex circuits, you would need to implement a more sophisticated
    #     evaluation function.

    #     Args:
    #         inputs: List of input polynomials, each represented as a list of coefficients
    #         modulus: The modulus for the polynomial arithmetic
    #         ring_dimension: The ring dimension for the polynomials

    #     Returns:
    #         List of output polynomials, each represented as a list of coefficients
    #     """
    #     # Check that we have the right number of inputs
    #     if len(inputs) != self.num_input:
    #         raise ValueError(f"Expected {self.num_input} inputs, got {len(inputs)}")

    #     # Initialize wires dictionary with input values
    #     wires = {}
    #     for i, input_poly in enumerate(inputs):
    #         wires[i] = input_poly

    #     # Process gates in order of their IDs
    #     for gate_id in sorted(self.gates.keys()):
    #         gate = self.gates[gate_id]

    #         if gate.gate_type == SerializablePolyGateType.INPUT:
    #             # Input gates are already handled
    #             continue

    #         elif gate.gate_type == SerializablePolyGateType.ADD:
    #             # Add two polynomials
    #             left = wires[gate.input_gates[0]]
    #             right = wires[gate.input_gates[1]]
    #             result = [(left[i] + right[i]) % modulus for i in range(ring_dimension)]
    #             wires[gate.gate_id] = result

    #         elif gate.gate_type == SerializablePolyGateType.SUB:
    #             # Subtract two polynomials
    #             left = wires[gate.input_gates[0]]
    #             right = wires[gate.input_gates[1]]
    #             result = [(left[i] - right[i]) % modulus for i in range(ring_dimension)]
    #             wires[gate.gate_id] = result

    #         elif gate.gate_type == SerializablePolyGateType.MUL:
    #             # Multiply two polynomials (simplified - just coefficient-wise multiplication)
    #             left = wires[gate.input_gates[0]]
    #             right = wires[gate.input_gates[1]]
    #             result = [(left[i] * right[i]) % modulus for i in range(ring_dimension)]
    #             wires[gate.gate_id] = result

    #         elif gate.gate_type == SerializablePolyGateType.SCALAR_MUL:
    #             # Multiply a polynomial by a scalar
    #             poly = wires[gate.input_gates[0]]
    #             # Use the coefficient array directly
    #             scalar_coeffs = gate.scalar_mul_data
    #             # Simplified scalar multiplication (just multiply by first coefficient)
    #             scalar = scalar_coeffs[0]
    #             result = [(poly[i] * scalar) % modulus for i in range(ring_dimension)]
    #             wires[gate.gate_id] = result

    #         elif gate.gate_type == SerializablePolyGateType.CALL:
    #             # Call a sub-circuit
    #             if gate.call_data["circuit_id"] not in self.sub_circuits:
    #                 raise ValueError(
    #                     f"Sub-circuit {gate.call_data['circuit_id']} not found"
    #                 )

    #             # Get the sub-circuit
    #             sub_circuit = self.sub_circuits[gate.call_data["circuit_id"]]

    #             # Prepare inputs for the sub-circuit
    #             sub_inputs = [wires[input_gate] for input_gate in gate.input_gates]

    #             # Evaluate the sub-circuit
    #             sub_outputs = sub_circuit.evaluate(sub_inputs, modulus, ring_dimension)

    #             # Store the output corresponding to this gate
    #             output_id = gate.call_data["output_id"]
    #             if output_id >= len(sub_outputs):
    #                 raise ValueError(
    #                     f"Output ID {output_id} out of range for sub-circuit outputs"
    #                 )
    #             wires[gate.gate_id] = sub_outputs[output_id]

    #     # Return the outputs
    #     return [wires[output_id] for output_id in self.output_ids]


# Example usage
if __name__ == "__main__":
    # Parameters for the polynomial
    ring_dimension = 8
    modulus = 257

    # Load a circuit from a JSON file
    circuit = SerializablePolyCircuit.load_from_file(
        "example_circuit.json", ring_dimension, modulus
    )

    # Example of how to access the circuit data
    print(f"Circuit has {len(circuit.gates)} gates")
    print(f"Circuit has {len(circuit.sub_circuits)} sub-circuits")
    print(f"Circuit has {circuit.num_input} inputs")
    print(f"Circuit has {len(circuit.output_ids)} outputs")

    # Example of how to access scalar multiplication data
    for gate_id, gate in circuit.gates.items():
        if gate.gate_type == SerializablePolyGateType.SCALAR_MUL:
            print(f"Gate {gate_id} is a scalar multiplication gate")
            print(
                f"Scalar data: {gate.scalar_mul_data[:10]}..."
            )  # Show first 10 coefficients
