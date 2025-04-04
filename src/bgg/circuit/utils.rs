use super::PolyCircuit;
use crate::bgg::circuit::{templates::ip_circuit, Evaluable};

/// Build a circuit that is a composition of two sub-circuits:
/// 1. The first sub-circuit is a public circuit that takes public inputs and outputs public outputs
/// 2. The second sub-circuit is an inner product circuit that takes private inputs and public
///    outputs and outputs their inner product
pub fn build_composite_circuit_from_public_and_ip<E: Evaluable>(
    public_circuit: PolyCircuit,
    num_priv_input: usize,
) -> PolyCircuit {
    let mut circuit = PolyCircuit::new();
    let num_pub_input = public_circuit.num_input();
    let num_input = num_pub_input + num_priv_input;
    let inputs = circuit.input(num_input);
    let pub_inputs = &inputs[0..num_pub_input];
    let priv_inputs = &inputs[num_pub_input..];

    // register the public circuit as sub-circuit
    let pub_circuit_id = circuit.register_sub_circuit(public_circuit);
    let pub_outputs = circuit.call_sub_circuit(pub_circuit_id, pub_inputs);

    // register the ip circuit as sub-circuit
    let ip_circuit = ip_circuit(priv_inputs.len(), pub_outputs.len());
    let ip_circuit_id = circuit.register_sub_circuit(ip_circuit);
    let ip_inputs = [priv_inputs.to_vec(), pub_outputs].concat();
    let ip_outputs = circuit.call_sub_circuit(ip_circuit_id, &ip_inputs);
    circuit.output(ip_outputs);
    circuit
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::{
        poly::{
            dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
            Poly,
        },
        utils::create_random_poly,
    };

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

        // Create a public circuit 
        let mut public_circuit = PolyCircuit::new();
        let pub_inputs = public_circuit.input(3);

        // Create 6 outputs (2 groups of 3)
        let out1 = public_circuit.add_gate(pub_inputs[0], pub_inputs[1]);
        let out2 = public_circuit.add_gate(pub_inputs[1], pub_inputs[2]);
        let out3 = public_circuit.add_gate(pub_inputs[0], pub_inputs[2]);
        let out4 = public_circuit.mul_gate(pub_inputs[0], pub_inputs[1]);
        let out5 = public_circuit.mul_gate(pub_inputs[1], pub_inputs[2]);
        let out6 = public_circuit.mul_gate(pub_inputs[0], pub_inputs[2]);

        public_circuit.output(vec![out1, out2, out3, out4, out5, out6]);

        // Create a copy of the public circuit for evaluation
        let public_circuit_for_eval = public_circuit.clone();

        // Evaluate the public circuit to get its outputs
        let pub_circuit_outputs =
            public_circuit_for_eval.eval(&params, &DCRTPoly::const_one(&params), &pub_polys);

        // Verify that the public circuit outputs are as expected
        let expected_out1 = pub_polys[0].clone() + pub_polys[1].clone(); // add_gate(pub_inputs[0], pub_inputs[1])
        let expected_out2 = pub_polys[1].clone() + pub_polys[2].clone(); // add_gate(pub_inputs[1], pub_inputs[2])
        let expected_out3 = pub_polys[0].clone() + pub_polys[2].clone(); // add_gate(pub_inputs[0], pub_inputs[2])
        let expected_out4 = &pub_polys[0] * &pub_polys[1]; // mul_gate(pub_inputs[0], pub_inputs[1])
        let expected_out5 = &pub_polys[1] * &pub_polys[2]; // mul_gate(pub_inputs[1], pub_inputs[2])
        let expected_out6 = &pub_polys[0] * &pub_polys[2]; // mul_gate(pub_inputs[0], pub_inputs[2])

        assert_eq!(pub_circuit_outputs.len(), 6);
        assert_eq!(pub_circuit_outputs[0], expected_out1);
        assert_eq!(pub_circuit_outputs[1], expected_out2);
        assert_eq!(pub_circuit_outputs[2], expected_out3);
        assert_eq!(pub_circuit_outputs[3], expected_out4);
        assert_eq!(pub_circuit_outputs[4], expected_out5);
        assert_eq!(pub_circuit_outputs[5], expected_out6);

        // Build the main circuit as a composition of the public circuit and the inner product
        // circuit
        let num_priv_input = 3;
        let main_circuit =
            build_composite_circuit_from_public_and_ip::<DCRTPoly>(public_circuit, num_priv_input);

        // Verify the circuit structure
        assert_eq!(main_circuit.num_input(), 6); // 3 private + 3 public inputs
        assert_eq!(main_circuit.num_output(), 2); // 2 inner products

        // Manually calculate the expected inner products
        // first inner product between the private inputs and the first 3 public outputs
        // second inner product between the private inputs and the last 3 public outputs
        let mut expected_ip1 = DCRTPoly::const_zero(&params);
        let mut expected_ip2 = DCRTPoly::const_zero(&params);

        for i in 0..num_priv_input {
            expected_ip1 += &pub_circuit_outputs[i] * &priv_polys[i];
            expected_ip2 += &pub_circuit_outputs[i + num_priv_input] * &priv_polys[i];
        }

        // Evaluate the main circuit
        let all_inputs = [pub_polys, priv_polys].concat();
        let result = main_circuit.eval(&params, &DCRTPoly::const_one(&params), &all_inputs);

        // Verify the results
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], expected_ip1);
        assert_eq!(result[1], expected_ip2);
    }
}
