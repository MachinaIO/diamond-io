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
    use num_bigint::BigUint;

    use super::*;
    use crate::{
        poly::{
            dcrt::{params::DCRTPolyParams, poly::DCRTPoly, DCRTPolyUniformSampler, FinRingElem},
            sampler::{DistType, PolyUniformSampler},
            Poly, PolyElem, PolyMatrix, PolyParams,
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

    #[test]
    fn test_fhe_eval_and_decrypt() {
        let params = DCRTPolyParams::default();
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let log_q = params.modulus_bits();

        // 1. Generate RLWE ciphertext (a, b) for input k
        let k = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
        let a = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let t = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let e = sampler_uniform.sample_uniform(&params, 1, 1, DistType::GaussDist { sigma: 0.0 });

        let modulus = params.modulus();
        let half_q = FinRingElem::half_q(&modulus.clone());
        let scale = DCRTPoly::from_const(&params, &half_q);
        let b = t.clone() * &a + &e - &(k.clone() * &scale);

        // decompose the polynomials a and b into bits
        let a_decomposed = a.get_column_matrix_decompose(0).get_column(0);
        let b_decomposed = b.get_column_matrix_decompose(0).get_column(0);

        assert!(a_decomposed.len() == b_decomposed.len());
        assert!(a_decomposed.len() == log_q);

        // 2. Create a public circuit with (a_decomposed, b_decomposed) hardcoded inside that takes
        // a public input constant polynomial x and returns 2*logq outputs: (a_decomposed[i] * x),
        // (b_decomposed[i] * x) for i = 0, 1, ..., logq - 1
        let mut public_circuit = PolyCircuit::new();
        let x_id = public_circuit.input(1)[0];

        // Convert a_decomposed and b_decomposed into a vector of bit vectors
        let a_bits_vecs: Vec<Vec<bool>> = a_decomposed
            .iter()
            .map(|poly| poly.coeffs().iter().map(|elem| elem.to_bit()).collect())
            .collect();

        let b_bits_vecs: Vec<Vec<bool>> = b_decomposed
            .iter()
            .map(|poly| poly.coeffs().iter().map(|elem| elem.to_bit()).collect())
            .collect();

        let mut public_circuit_outputs = Vec::new();

        for i in 0..a_bits_vecs.len() {
            let a_ith_bits_const = public_circuit.const_bit_poly(&a_bits_vecs[i]);
            let b_ith_bits_const = public_circuit.const_bit_poly(&b_bits_vecs[i]);
            let a_ith_bits_and_x = public_circuit.mul_gate(a_ith_bits_const, x_id);
            let b_ith_bits_and_x = public_circuit.mul_gate(b_ith_bits_const, x_id);
            public_circuit_outputs.push(a_ith_bits_and_x);
            public_circuit_outputs.push(b_ith_bits_and_x);
        }

        public_circuit.output(public_circuit_outputs);

        // Create a copy of the public circuit for evaluation
        let public_circuit_copy = public_circuit.clone();

        // Evaluate the copied public circuit
        let x = DCRTPoly::const_one(&params);
        let public_output =
            public_circuit_copy.eval(&params, &DCRTPoly::const_one(&params), &[x.clone()]);

        // Check output length
        assert_eq!(public_output.len(), 2 * log_q);

        // Compute expected output manually
        let mut expected_output = Vec::new();
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_decomposed[i].clone();
            let b_ith_bits = b_decomposed[i].clone();
            expected_output.push(a_ith_bits * x.clone());
            expected_output.push(b_ith_bits * x.clone());
        }

        // Compare each output
        for (actual, expected) in public_output.iter().zip(expected_output.iter()) {
            assert_eq!(actual, expected);
        }

        // 3. Build the main circuit as a composition of the public circuit and the inner product
        // circuit such that the inner product circuit performs the inner product of the public
        // circuit outputs with the private inputs [-t, 1]
        let main_circuit =
            build_composite_circuit_from_public_and_ip::<DCRTPoly>(public_circuit, 2);

        // Verify the circuit structure
        assert_eq!(main_circuit.num_input(), 3); // 1 public inputs + 2 private inputs
        assert_eq!(main_circuit.num_output(), log_q); // logq inner products

        // Evaluate the main circuit
        let all_inputs = [x.clone(), -(t.entry(0, 0)), DCRTPoly::const_one(&params)];
        let decomposed_result =
            main_circuit.eval(&params, &DCRTPoly::const_one(&params), &all_inputs);

        // Verify the results
        assert_eq!(decomposed_result.len(), log_q);

        for i in 0..log_q {
            let expected_inner_product = (a_decomposed[i].clone() * x.clone()) * -(t.entry(0, 0)) +
                (b_decomposed[i].clone() * x.clone());
            assert_eq!(decomposed_result[i], expected_inner_product);
        }

        let recomposed_result = DCRTPoly::from_decomposed(&params, &decomposed_result);

        let prod = x.clone() * k.entry(0, 0);
        let prod_bits = prod.coeffs().iter().map(|elem| elem.to_bit()).collect::<Vec<_>>();

        for (i, coeff) in recomposed_result.coeffs().iter().enumerate() {
            let actual = {
                if coeff.value() != &BigUint::from(0u8) {
                    true
                } else {
                    false
                }
            };
            let expected = prod_bits[i];
            assert_eq!(actual, expected);
        }
    }
}
