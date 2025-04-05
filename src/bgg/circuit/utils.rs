use super::PolyCircuit;
use crate::bgg::circuit::Evaluable;

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
    // - outputs: b'_decomposed[i] - a'_decomposed[i] * t - (for i = 0, 1, ..., logq - 1)
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
        let b_i_minus_a_i_times_t = dec_circuit.sub_gate(b_i, a_i_times_t);
        lhs.push(b_i_minus_a_i_times_t);
    }
    dec_circuit.output(lhs.clone());

    // register the decryption circuit as sub-circuit
    let dec_circuit_id = circuit.register_sub_circuit(dec_circuit);
    let dec_inputs = [pub_circuit_outputs.clone(), vec![t_in]].concat();
    let dec_outputs = circuit.call_sub_circuit(dec_circuit_id, &dec_inputs);

    circuit.output(dec_outputs);
    circuit
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{params::DCRTPolyParams, poly::DCRTPoly, DCRTPolyUniformSampler, FinRingElem},
        sampler::DistType,
        Poly, PolyElem, PolyParams,
    };

    #[test]
    fn test_fhe_eval_and_decrypt() {
        let params = DCRTPolyParams::default();
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let log_q = params.modulus_bits();

        // 1. Generate RLWE ciphertext (a, b) for input k
        // b = a * t - k * q/2 + e
        let k = sampler_uniform.sample_poly(&params, &DistType::BitDist);
        let t = sampler_uniform.sample_poly(&params, &DistType::FinRingDist);
        let a = sampler_uniform.sample_poly(&params, &DistType::FinRingDist);
        let e = sampler_uniform.sample_poly(&params, &DistType::GaussDist { sigma: 0.0 });

        let modulus = params.modulus();
        let half_q = FinRingElem::half_q(&modulus.clone());
        let scale = DCRTPoly::from_const(&params, &half_q);
        let b = a.clone() * t.clone() - &(k.clone() * &scale) + &e;

        // decompose the polynomials a and b into bits
        let a_bits = a.decompose(&params);
        let b_bits = b.decompose(&params);

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
        let all_inputs = [x.clone(), t.clone()];
        let output = main_circuit.eval(&params, &DCRTPoly::const_one(&params), &all_inputs);

        // Verify the results
        assert_eq!(output.len(), log_q);

        // Compute expected output manually
        let mut expected_output = Vec::new();
        for i in 0..a_bits_vecs.len() {
            let a_ith_bits = a_bits[i].clone();
            let b_ith_bits = b_bits[i].clone();
            expected_output.push((b_ith_bits * x.clone()) - (a_ith_bits * x.clone()) * t.clone());
            assert_eq!(output[i], expected_output[i]);
        }

        // Recompose the output
        let output_recomposed = DCRTPoly::from_decomposed(&params, &output);

        // Compute decision threshold values
        let modulus = params.modulus();
        let quarter_q = modulus.as_ref() >> 2; // q/4
        let three_quarter_q = &quarter_q * 3u32; // 3q/4

        // Decode plaintext directly into boolean vector
        let recovered_bits: Vec<bool> = output_recomposed
            .coeffs()
            .iter()
            .map(|coeff| coeff.value())
            .map(|coeff| coeff > &quarter_q && coeff <= &three_quarter_q)
            .collect();

        // Verify correctness
        assert_eq!(recovered_bits.len(), params.ring_dimension() as usize);
        assert_eq!(recovered_bits, (k * x).to_bool_vec());
    }
}
