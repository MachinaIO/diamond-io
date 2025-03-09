use super::{Evaluable, PolyCircuit};
use crate::poly::*;

/// Build a circuit that takes as input a public circuit and a private input, and outputs inner products between the `num_priv_input`-sized output of the public circuit and the private input.
pub fn build_ip_priv_and_pub_circuit_outputs<P: Poly, E: Evaluable<P>>(
    public_circuit: PolyCircuit<P>,
    num_priv_input: usize,
) -> PolyCircuit<P> {
    let num_pub_input = public_circuit.num_input();
    debug_assert_eq!(public_circuit.num_output() % num_priv_input, 0);
    let num_ip_outputs = public_circuit.num_output() / num_priv_input;
    let num_input = num_pub_input + num_priv_input;
    let mut circuit = PolyCircuit::<P>::new();
    let inputs = circuit.input(num_input);
    let priv_inputs = &inputs[0..num_priv_input];
    let pub_inputs = &inputs[num_priv_input..];
    let circuit_id = circuit.register_sub_circuit(public_circuit);
    let pub_outputs = circuit.call_sub_circuit(circuit_id, pub_inputs);
    let mut ip_outputs = Vec::new();
    for out_idx in 0..num_ip_outputs {
        let mut ip_output: Option<usize> = None;
        for (priv_input, pub_output) in priv_inputs
            .iter()
            .zip(pub_outputs[(num_priv_input * out_idx)..(num_priv_input * (out_idx + 1))].iter())
        {
            let mul = circuit.mul_gate(*pub_output, *priv_input);
            if let Some(prev_ip) = ip_output {
                ip_output = Some(circuit.add_gate(prev_ip, mul));
            } else {
                ip_output = Some(mul);
            }
        }
        ip_outputs.push(ip_output.unwrap());
    }
    circuit.output(ip_outputs);
    circuit
}

/// Build a circuit that takes as input `num_bits` bits and outputs the integer represented by the bits.
pub fn build_bits_to_int_circuit<P: Poly, E: Evaluable<P>>(
    params: &P::Params,
    num_bits: usize,
) -> PolyCircuit<P> {
    let mut circuit = PolyCircuit::<P>::new();
    let inputs = circuit.input(num_bits);
    let mut output = None;
    for idx in 0..num_bits {
        let scalar_poly = P::const_power_of_two(params, idx);
        let scalar_muled = circuit.scalar_mul_gate(inputs[idx], scalar_poly);
        if let Some(prev_output) = output {
            output = Some(circuit.add_gate(prev_output, scalar_muled));
        } else {
            output = Some(scalar_muled);
        }
    }
    circuit.output(vec![output.unwrap()]);
    circuit
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;

    use super::*;
    use crate::poly::dcrt::FinRingElem;
    use crate::poly::dcrt::{
        params::DCRTPolyParams, poly::DCRTPoly, sampler::uniform::DCRTPolyUniformSampler,
    };
    use crate::poly::sampler::DistType;

    // Helper function to create a random polynomial using UniformSampler
    fn create_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
        let sampler = DCRTPolyUniformSampler::new();
        sampler.sample_poly(params, &DistType::FinRingDist)
    }

    // Helper function to create a bit polynomial (0 or 1)
    fn create_bit_poly(params: &DCRTPolyParams, bit: bool) -> DCRTPoly {
        if bit {
            DCRTPoly::const_one(params)
        } else {
            DCRTPoly::const_zero(params)
        }
    }

    #[test]
    fn test_build_ip_priv_and_pub_circuit_outputs() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Create input polynomials
        let priv_polys = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        let pub_polys = vec![
            create_random_poly(&params),
            create_random_poly(&params),
            create_random_poly(&params),
        ];

        // Create a simple public circuit that adds its inputs
        let mut public_circuit = PolyCircuit::<DCRTPoly>::new();
        let pub_inputs = public_circuit.input(3);

        // Create 6 outputs (2 groups of 3)
        let out1 = public_circuit.add_gate(pub_inputs[0], pub_inputs[1]);
        let out2 = public_circuit.add_gate(pub_inputs[1], pub_inputs[2]);
        let out3 = public_circuit.add_gate(pub_inputs[0], pub_inputs[2]);
        let out4 = public_circuit.mul_gate(pub_inputs[0], pub_inputs[1]);
        let out5 = public_circuit.mul_gate(pub_inputs[1], pub_inputs[2]);
        let out6 = public_circuit.mul_gate(pub_inputs[0], pub_inputs[2]);

        public_circuit.output(vec![out1, out2, out3, out4, out5, out6]);

        // Number of private inputs (we'll have 2 inner products)
        let num_priv_input = 3;

        // Create a copy of the public circuit for evaluation
        let mut public_circuit_for_eval = PolyCircuit::<DCRTPoly>::new();
        let pub_inputs_eval = public_circuit_for_eval.input(3);

        let out1_eval = public_circuit_for_eval.add_gate(pub_inputs_eval[0], pub_inputs_eval[1]);
        let out2_eval = public_circuit_for_eval.add_gate(pub_inputs_eval[1], pub_inputs_eval[2]);
        let out3_eval = public_circuit_for_eval.add_gate(pub_inputs_eval[0], pub_inputs_eval[2]);
        let out4_eval = public_circuit_for_eval.mul_gate(pub_inputs_eval[0], pub_inputs_eval[1]);
        let out5_eval = public_circuit_for_eval.mul_gate(pub_inputs_eval[1], pub_inputs_eval[2]);
        let out6_eval = public_circuit_for_eval.mul_gate(pub_inputs_eval[0], pub_inputs_eval[2]);

        public_circuit_for_eval
            .output(vec![out1_eval, out2_eval, out3_eval, out4_eval, out5_eval, out6_eval]);

        // Evaluate the public circuit to get its outputs
        let pub_circuit_outputs = public_circuit_for_eval.eval_poly_circuit(
            &params,
            DCRTPoly::const_one(&params),
            &pub_polys,
        );

        // Verify that the public circuit outputs are as expected
        let expected_out1 = pub_polys[0].clone() + pub_polys[1].clone(); // add_gate(pub_inputs[0], pub_inputs[1])
        let expected_out2 = pub_polys[1].clone() + pub_polys[2].clone(); // add_gate(pub_inputs[1], pub_inputs[2])
        let expected_out3 = pub_polys[0].clone() + pub_polys[2].clone(); // add_gate(pub_inputs[0], pub_inputs[2])
        let expected_out4 = pub_polys[0].clone() * pub_polys[1].clone(); // mul_gate(pub_inputs[0], pub_inputs[1])
        let expected_out5 = pub_polys[1].clone() * pub_polys[2].clone(); // mul_gate(pub_inputs[1], pub_inputs[2])
        let expected_out6 = pub_polys[0].clone() * pub_polys[2].clone(); // mul_gate(pub_inputs[0], pub_inputs[2])

        assert_eq!(pub_circuit_outputs.len(), 6);
        assert_eq!(pub_circuit_outputs[0], expected_out1);
        assert_eq!(pub_circuit_outputs[1], expected_out2);
        assert_eq!(pub_circuit_outputs[2], expected_out3);
        assert_eq!(pub_circuit_outputs[3], expected_out4);
        assert_eq!(pub_circuit_outputs[4], expected_out5);
        assert_eq!(pub_circuit_outputs[5], expected_out6);

        // Build the inner product circuit
        let ip_circuit = build_ip_priv_and_pub_circuit_outputs::<DCRTPoly, DCRTPoly>(
            public_circuit,
            num_priv_input,
        );

        // Verify the circuit structure
        assert_eq!(ip_circuit.num_input(), 6); // 3 private + 3 public inputs
        assert_eq!(ip_circuit.num_output(), 2); // 2 inner products

        // Manually calculate the expected inner products
        let mut expected_ip1 = DCRTPoly::const_zero(&params);
        let mut expected_ip2 = DCRTPoly::const_zero(&params);

        for i in 0..num_priv_input {
            expected_ip1 = expected_ip1 + (pub_circuit_outputs[i].clone() * priv_polys[i].clone());
            expected_ip2 = expected_ip2
                + (pub_circuit_outputs[i + num_priv_input].clone() * priv_polys[i].clone());
        }

        // Evaluate the inner product circuit
        let all_inputs = [priv_polys, pub_polys].concat();
        let result =
            ip_circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &all_inputs);

        // Verify the results
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], expected_ip1);
        assert_eq!(result[1], expected_ip2);
    }

    #[test]
    fn test_bits_to_int_circuit() {
        // Create parameters for testing
        let params = DCRTPolyParams::default();

        // Test with 4 bits
        let num_bits = 4;
        let circuit = build_bits_to_int_circuit::<DCRTPoly, DCRTPoly>(&params, num_bits);

        // Verify the circuit structure
        assert_eq!(circuit.num_input(), num_bits);
        assert_eq!(circuit.num_output(), 1);

        // Test case 1: 0101 (binary) = 5 (decimal)
        let bits1 = vec![
            create_bit_poly(&params, true),  // 1 (LSB)
            create_bit_poly(&params, false), // 0
            create_bit_poly(&params, true),  // 1
            create_bit_poly(&params, false), // 0 (MSB)
        ];

        let result1 = circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &bits1);

        // Expected: 1*2^0 + 0*2^1 + 1*2^2 + 0*2^3 = 5
        let expected1 = DCRTPoly::from_const(&params, &FinRingElem::new(5, params.modulus()));

        assert_eq!(result1.len(), 1);
        assert_eq!(result1[0], expected1);

        // Test case 2: 1111 (binary) = 15 (decimal)
        let bits2 = vec![
            create_bit_poly(&params, true), // 1 (LSB)
            create_bit_poly(&params, true), // 1
            create_bit_poly(&params, true), // 1
            create_bit_poly(&params, true), // 1 (MSB)
        ];

        let result2 = circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &bits2);

        // Expected: 1*2^0 + 1*2^1 + 1*2^2 + 1*2^3 = 15
        let expected2 = DCRTPoly::from_const(&params, &FinRingElem::new(15, params.modulus()));

        assert_eq!(result2.len(), 1);
        assert_eq!(result2[0], expected2);

        // Test case 3: 0000 (binary) = 0 (decimal)
        let bits3 = vec![
            create_bit_poly(&params, false), // 0 (LSB)
            create_bit_poly(&params, false), // 0
            create_bit_poly(&params, false), // 0
            create_bit_poly(&params, false), // 0 (MSB)
        ];

        let result3 = circuit.eval_poly_circuit(&params, DCRTPoly::const_one(&params), &bits3);

        // Expected: 0
        let expected3 = DCRTPoly::const_zero(&params);

        assert_eq!(result3.len(), 1);
        assert_eq!(result3[0], expected3);
    }
}
