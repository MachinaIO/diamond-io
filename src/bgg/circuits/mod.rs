pub mod circuit;
pub mod gate;
pub mod plain_eval;

pub use circuit::*;
pub use gate::*;
pub use plain_eval::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{params::DCRTPolyParams, poly::DCRTPoly, sampler::uniform::DCRTPolyUniformSampler},
        sampler::{DistType, PolyUniformSampler},
    };

    // Helper function to create a random polynomial using UniformSampler
    fn create_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
        let sampler = DCRTPolyUniformSampler::new();
        sampler.sample_poly(params, &DistType::FinRingDist)
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
        let result = eval_poly_circuit(circuit, &[poly1.clone(), poly2.clone()]);

        // Expected result: poly1 + poly2
        let expected = poly1 + poly2;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
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
        let result = eval_poly_circuit(circuit, &[poly1.clone(), poly2.clone()]);

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
        let result = eval_poly_circuit(circuit, &[poly.clone()]);

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
        let result = eval_poly_circuit(circuit, &[poly1.clone(), poly2.clone()]);

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
        let result = eval_poly_circuit(circuit, &[poly1.clone(), poly2.clone(), poly3.clone()]);

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
        let result = eval_poly_circuit(circuit, &[poly1.clone(), poly2.clone()]);

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
        let result = eval_poly_circuit(
            circuit,
            &[poly1.clone(), poly2.clone(), poly3.clone(), poly4.clone()],
        );

        // Expected result: (((poly1 + poly2) * (poly3 * poly4)) + (poly1 - poly3)) * scalar
        let expected =
            (((poly1.clone() + poly2) * (poly3.clone() * poly4)) + (poly1 - poly3)) * scalar;

        // Verify the result
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], expected);
    }
}
