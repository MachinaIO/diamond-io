use super::PolyCircuit;
use crate::bgg::circuit::Evaluable;

/// Build a circuit that is a composition of two sub-circuits:
/// 1. A public circuit that it is assumed to return ciphertexts where each ciphertext is bit
///    decomposed (a_decomposed[0], ... a_decomposed[logq-1] , b_decomposed[0], ...
///    b_decomposed[logq-1])
/// 2. An FHE decryption circuit that takes each ciphertext and the RLWE secret key t as inputs and
///    returns the bit decomposed plaintext for each cipheretxt
pub fn build_composite_circuit_from_public_and_fhe_dec<E: Evaluable>(
    public_circuit: PolyCircuit,
    log_q: usize,
) -> PolyCircuit {
    let num_input_pub_circuit = public_circuit.num_input();
    let num_output_pub_circuit = public_circuit.num_output();
    assert_eq!(num_output_pub_circuit % (2 * log_q), 0);

    let mut circuit = PolyCircuit::new();
    let inputs = circuit.input(num_input_pub_circuit + 1); // the extra input is the secret key t
    let t_in = inputs[num_input_pub_circuit];
    let pub_circuit_inputs = &inputs[0..num_input_pub_circuit];

    // register the public circuit as sub-circuit
    let pub_circuit_id = circuit.register_sub_circuit(public_circuit);
    let pub_circuit_outputs = circuit.call_sub_circuit(pub_circuit_id, pub_circuit_inputs);

    let num_ciphertexts = pub_circuit_outputs.len() / (2 * log_q);

    // build the FHE decryption circuit:
    // - inputs: a_decomposed[0], ... a_decomposed[logq-1] , b_decomposed[0], ...
    //   b_decomposed[logq-1] for each ciphertext, t ((2 * logq * num_ciphertexts + 1 inputs))
    // - outputs: b'_decomposed[i] - a'_decomposed[i] * t - (for i = 0, 1, ..., logq - 1) for each
    //   ciphertext
    // (logq * num_ciphertexts outputs)
    let mut dec_circuit = PolyCircuit::new();
    let dec_circuit_inputs = dec_circuit.input(num_output_pub_circuit + 1);
    let mut output = Vec::with_capacity(log_q * num_ciphertexts);
    let t = dec_circuit_inputs[num_output_pub_circuit];
    for ct_idx in 0..num_ciphertexts {
        let offset = ct_idx * log_q;
        let a_start = offset;
        let b_start = log_q * num_ciphertexts + offset;
        for bit_idx in 0..log_q {
            let a_i = dec_circuit_inputs[a_start + bit_idx];
            let b_i = dec_circuit_inputs[b_start + bit_idx];
            let a_i_times_t = dec_circuit.mul_gate(a_i, t);
            let b_i_minus_a_i_times_t = dec_circuit.sub_gate(b_i, a_i_times_t);
            output.push(b_i_minus_a_i_times_t);
        }
    }
    dec_circuit.output(output);

    // register the decryption circuit as sub-circuit
    let dec_circuit_id = circuit.register_sub_circuit(dec_circuit);
    let dec_inputs = [pub_circuit_outputs, vec![t_in]].concat();
    let dec_outputs = circuit.call_sub_circuit(dec_circuit_id, &dec_inputs);

    circuit.output(dec_outputs);
    circuit
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::{
        io::utils::encrypt_rlwe,
        poly::{
            dcrt::{params::DCRTPolyParams, poly::DCRTPoly, DCRTPolyUniformSampler},
            sampler::DistType,
            Poly, PolyParams,
        },
    };

    #[test]
    fn test_fhe_eval_and_decrypt() {
        let params = DCRTPolyParams::default();
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let log_q = params.modulus_bits();
        let sigma = 3.0;

        // 1. Generate RLWE ciphertext (a, b) for input k
        // b = a * t - k * q/2 + e
        let k = sampler_uniform.sample_poly(&params, &DistType::BitDist);
        let (a, b, t) = encrypt_rlwe(&params, &sampler_uniform, sigma, &k);

        // decompose the polynomials a and b into bits
        let a_bits = a.decompose(&params);
        let b_bits = b.decompose(&params);

        assert!(a_bits.len() == b_bits.len());
        assert!(a_bits.len() == log_q);

        // 2. Create a public circuit
        // Input: a_bits[0], ..., a_bits[logq - 1], b_bits[0], ..., b_bits[logq - 1], x, y
        // Output: a_bits[0] AND x, ..., a_bits[logq - 1] AND x, a_bits[0] AND y, ..., a_bits[logq -
        // 1] AND y, b_bits[0] AND x, ..., b_bits[logq - 1] AND x, b_bits[0] AND y, ..., b_bits[logq
        // - 1] AND y
        // a_bits and b_bits are hardcoded inside the circuit

        let a_bits_vecs = a_bits.iter().map(|poly| poly.to_bool_vec()).collect::<Vec<_>>();
        let b_bits_vecs = b_bits.iter().map(|poly| poly.to_bool_vec()).collect::<Vec<_>>();
        let mut public_circuit = PolyCircuit::new();
        let pub_circuit_inputs = public_circuit.input(2);
        let x = pub_circuit_inputs[0];
        let y = pub_circuit_inputs[1];

        let mut public_circuit_outputs = Vec::new();

        let a_bits_consts =
            a_bits_vecs.iter().map(|poly| public_circuit.const_bit_poly(poly)).collect::<Vec<_>>();

        let b_bits_consts =
            b_bits_vecs.iter().map(|poly| public_circuit.const_bit_poly(poly)).collect::<Vec<_>>();

        for i in 0..a_bits_vecs.len() {
            let a_ith_bits_and_x = public_circuit.and_gate(a_bits_consts[i], x);
            public_circuit_outputs.push(a_ith_bits_and_x);
        }

        for i in 0..a_bits_vecs.len() {
            let a_ith_bits_and_y = public_circuit.and_gate(a_bits_consts[i], y);
            public_circuit_outputs.push(a_ith_bits_and_y);
        }

        for i in 0..b_bits_vecs.len() {
            let b_ith_bits_and_x = public_circuit.and_gate(b_bits_consts[i], x);
            public_circuit_outputs.push(b_ith_bits_and_x);
        }

        for i in 0..b_bits_vecs.len() {
            let b_ith_bits_and_y = public_circuit.and_gate(b_bits_consts[i], y);
            public_circuit_outputs.push(b_ith_bits_and_y);
        }

        public_circuit.output(public_circuit_outputs);

        // Create a copy of the public circuit for evaluation
        let public_circuit_copy = public_circuit.clone();

        // Evaluate the copied public circuit
        let x = DCRTPoly::const_one(&params);
        let y = DCRTPoly::const_zero(&params);
        let public_output = public_circuit_copy.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[x.clone(), y.clone()],
        );

        // Check output length
        assert_eq!(public_circuit_copy.num_input(), 2);
        assert_eq!(public_output.len(), (2 * log_q) * 2);

        // Compute expected output manually and verify
        let mut expected_output = Vec::new();
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            expected_output.push(a_ith_bits * x.clone());
        }
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            expected_output.push(a_ith_bits * y.clone());
        }
        for i in 0..b_bits_vecs.len() {
            let b_ith_bits = b_bits[i].clone();
            expected_output.push(b_ith_bits * x.clone());
        }
        for i in 0..b_bits_vecs.len() {
            let b_ith_bits = b_bits[i].clone();
            expected_output.push(b_ith_bits * y.clone());
        }

        for (actual, expected) in public_output.iter().zip(expected_output.iter()) {
            assert_eq!(actual, expected);
        }

        // 3. Build the main circuit as a composition of the public circuit and the FHE decryption
        //    circuit
        let main_circuit =
            build_composite_circuit_from_public_and_fhe_dec::<DCRTPoly>(public_circuit, log_q);

        // Verify the circuit structure
        assert_eq!(main_circuit.num_input(), 3); // 2 public inputs (x, y) and 1 private input (t)
        assert_eq!(main_circuit.num_output(), 2 * log_q); // 2 * log_q outputs

        // Evaluate the main circuit
        let all_inputs = [x.clone(), y.clone(), t.clone()];
        let output = main_circuit.eval(&params, &DCRTPoly::const_one(&params), &all_inputs);

        // Verify the results
        assert_eq!(output.len(), 2 * log_q);

        // Compute expected output manually
        let mut expected_output = Vec::new();
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            let b_ith_bits = b_bits[i].clone();
            expected_output.push((b_ith_bits * x.clone()) - (a_ith_bits * x.clone()) * t.clone());
        }
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            let b_ith_bits = b_bits[i].clone();
            expected_output.push((b_ith_bits * y.clone()) - (a_ith_bits * y.clone()) * t.clone());
        }

        // Verify the output
        for (actual, expected) in output.iter().zip(expected_output.iter()) {
            assert_eq!(actual, expected);
        }
        // 5. Recompose the outputs
        let output_one = output[..log_q].to_vec();
        let output_two = output[log_q..].to_vec();

        let output_one_recomposed = DCRTPoly::from_decomposed(&params, &output_one);
        let output_two_recomposed = DCRTPoly::from_decomposed(&params, &output_two);

        // recover the bits
        let recovered_bits_one = output_one_recomposed.extract_bits_with_threshold(&params);
        let recovered_bits_two = output_two_recomposed.extract_bits_with_threshold(&params);

        // Verify correctness
        assert_eq!(recovered_bits_one.len(), params.ring_dimension() as usize);
        assert_eq!(recovered_bits_two.len(), params.ring_dimension() as usize);
        assert_eq!(recovered_bits_one, (k.clone() * x).to_bool_vec());
        assert_eq!(recovered_bits_two, (k * y).to_bool_vec());
    }
}
