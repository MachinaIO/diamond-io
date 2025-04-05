use super::PolyCircuit;
use crate::bgg::circuit::{templates::ip_circuit, Evaluable};

/// Build a circuit that is a composition of two sub-circuits:
/// 1. A public circuit that it is assumed to return a ciphertext in its bit decomposed form
///    (a_decomposed[0], ... a_decomposed[logq-1] , b_decomposed[0], ... b_decomposed[logq-1])
/// 2. An FHE decryption circuit that takes the ciphertext and the RLWE secret key t as inputs and
///    returns the bit decomposed plaintext
pub fn build_composite_circuit_from_public_and_fhe_dec<E: Evaluable>(
    public_circuit: PolyCircuit,
) -> PolyCircuit {
    let num_input_pub_circuit = public_circuit.num_input();

    let mut circuit = PolyCircuit::new();
    let inputs = circuit.input(num_input_pub_circuit + 1); // the extra input is the secret key t
    let t_in = inputs[num_input_pub_circuit];
    let pub_circuit_inputs = &inputs[0..num_input_pub_circuit];

    // register the public circuit as sub-circuit
    let pub_circuit_id = circuit.register_sub_circuit(public_circuit);
    let pub_circuit_outputs = circuit.call_sub_circuit(pub_circuit_id, pub_circuit_inputs);

    // build the FHE decryption circuit:
    // - inputs: a_decomposed[0], ... a_decomposed[logq-1] , b_decomposed[0], ...
    //   b_decomposed[logq-1], t (2 * logq + 1 inputs)
    // - outputs: a'_decomposed[i] * t - b'_decomposed[i] (for i = 0, 1, ..., logq - 1)
    // (logq outputs)
    let mut dec_circuit = PolyCircuit::new();
    let dec_circuit_inputs = dec_circuit.input(pub_circuit_outputs.len() + 1);
    let num_output_dec_circuit = pub_circuit_outputs.len() / 2;
    let mut lhs = Vec::with_capacity(num_output_dec_circuit);
    let rhs = dec_circuit_inputs[2 * num_output_dec_circuit];
    for i in 0..num_output_dec_circuit {
        let a_i = dec_circuit_inputs[i];
        let b_i = dec_circuit_inputs[i + num_output_dec_circuit];
        let a_i_times_t = dec_circuit.mul_gate(a_i, rhs);
        let a_i_t_minus_b_i = dec_circuit.sub_gate(a_i_times_t, b_i);
        lhs.push(a_i_t_minus_b_i);
    }
    dec_circuit.output(lhs.clone());

    // register the decryption circuit as sub-circuit
    let dec_circuit_id = circuit.register_sub_circuit(dec_circuit);
    let dec_inputs = [pub_circuit_outputs.clone(), vec![t_in]].concat();
    let dec_outputs = circuit.call_sub_circuit(dec_circuit_id, &dec_inputs);

    circuit.output(dec_outputs);
    circuit
}

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
        // b = a * t - k * q/2 + e
        let k = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
        let t = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let a = sampler_uniform.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let e = sampler_uniform.sample_uniform(&params, 1, 1, DistType::GaussDist { sigma: 0.0 });

        let modulus = params.modulus();
        let half_q = FinRingElem::half_q(&modulus.clone());
        let scale = DCRTPoly::from_const(&params, &half_q);
        let b = a.clone() * t.clone() - &(k.clone() * &scale) + &e;

        // decompose the polynomials a and b into bits
        let a_bits = a.get_column_matrix_decompose(0).get_column(0);
        let b_bits = b.get_column_matrix_decompose(0).get_column(0);

        assert!(a_bits.len() == b_bits.len());
        assert!(a_bits.len() == log_q);

        // 2. Create a public circuit
        // Input: a_bits[0], ..., a_bits[logq - 1], b_bits[0], ..., b_bits[logq - 1], x
        // Output: a_bits[0] AND x, ..., a_bits[logq - 1] AND x, b_bits[0] AND x, ..., b_bits[logq -
        // 1] AND x
        // a_bits and b_bits are hardcoded inside the circuit

        let a_bits_vecs = a_bits.iter().map(|poly| poly.to_bool_vec()).collect::<Vec<_>>();
        let b_bits_vecs = b_bits.iter().map(|poly| poly.to_bool_vec()).collect::<Vec<_>>();
        let mut public_circuit = PolyCircuit::new();
        let x_id = public_circuit.input(1)[0];

        let mut public_circuit_outputs = Vec::new();

        for i in 0..a_bits_vecs.len() {
            let a_ith_bits_const = public_circuit.const_bit_poly(&a_bits_vecs[i]);
            let a_ith_bits_and_x = public_circuit.and_gate(a_ith_bits_const, x_id);
            public_circuit_outputs.push(a_ith_bits_and_x);
        }

        for i in 0..b_bits_vecs.len() {
            let b_ith_bits_const = public_circuit.const_bit_poly(&b_bits_vecs[i]);
            let b_ith_bits_and_x = public_circuit.and_gate(b_ith_bits_const, x_id);
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
        assert_eq!(public_circuit_copy.num_input(), 1);
        assert_eq!(public_output.len(), 2 * log_q);

        // Compute expected output manually
        let mut expected_output = Vec::new();
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            expected_output.push(a_ith_bits * x.clone());
        }
        for i in 0..b_bits_vecs.len() {
            let b_ith_bits = b_bits[i].clone();
            expected_output.push(b_ith_bits * x.clone());
        }

        // Compare each output
        for (actual, expected) in public_output.iter().zip(expected_output.iter()) {
            assert_eq!(actual, expected);
        }

        // 3. Build the main circuit as a composition of the public circuit and the FHE decryption
        //    circuit
        let main_circuit =
            build_composite_circuit_from_public_and_fhe_dec::<DCRTPoly>(public_circuit);

        // Verify the circuit structure
        assert_eq!(main_circuit.num_input(), 2); // 1 public input (x) and 1 private input (t)
        assert_eq!(main_circuit.num_output(), log_q); // log_q outputs

        // Evaluate the main circuit
        let all_inputs = [x.clone(), t.entry(0, 0)];
        let output = main_circuit.eval(&params, &DCRTPoly::const_one(&params), &all_inputs);

        // Verify the results
        assert_eq!(output.len(), log_q);

        // Compute expected output manually
        let mut expected_output = Vec::new();
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            let b_ith_bits = b_bits[i].clone();
            expected_output
                .push((a_ith_bits * x.clone()) * t.entry(0, 0).clone() - (b_ith_bits * x.clone()));
            assert_eq!(output[i], expected_output[i]);
        }

        // Recompose the output
        let output_recomposed = DCRTPoly::from_decomposed(&params, &output);

        // decrypt the result
        let plaintext_recovered = (output_recomposed).extract_highest_bits();
        assert_eq!(plaintext_recovered.len(), params.ring_dimension() as usize);
        assert_eq!(plaintext_recovered, (k.entry(0, 0) * x.clone()).to_bool_vec());
    }
}
