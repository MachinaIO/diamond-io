use crate::bgg::circuit::Evaluable;

use super::PolyCircuit;

/// Build a circuit that is a composition of two sub-circuits:
/// 1. A public circuit that it is assumed to return one or more ciphertexts where each ciphertext
///    is bit decomposed -> BITS(ct_0), BITS(ct_1), ... where BITS(ct_i) = [a_bit_0, b_bit_0,
///    a_bit_1, b_bit_1, ...]
/// 2. An FHE decryption circuit that takes each ciphertext and the RLWE secret key -t_bar as inputs
///    and returns the bit decomposed plaintext for each cipheretxt
pub fn build_composite_circuit_from_public_and_fhe_dec<E: Evaluable>(
    public_circuit: PolyCircuit,
    log_q: usize,
) -> PolyCircuit {
    let num_pub_circuit_input = public_circuit.num_input();
    let num_pub_circuit_output = public_circuit.num_output();
    debug_assert_eq!(num_pub_circuit_output % (2 * log_q), 0);

    let num_output = public_circuit.num_output() / 2;
    let num_input = num_pub_circuit_input + 1;
    let mut circuit = PolyCircuit::new();
    let inputs = circuit.input(num_input);
    let pub_circuit_inputs = &inputs[0..num_pub_circuit_input];
    let minus_t_bar = &inputs[num_pub_circuit_input];
    let circuit_id = circuit.register_sub_circuit(public_circuit);
    let pub_circuit_outputs = circuit.call_sub_circuit(circuit_id, pub_circuit_inputs);
    let mut outputs = Vec::with_capacity(num_output);

    // Process each pair (a_bit, b_bit) from the public circuit outputs
    // result = a_bit * -t_bar + b_bit
    for i in 0..num_output {
        let a_bit = pub_circuit_outputs[i * 2];
        let b_bit = pub_circuit_outputs[i * 2 + 1];
        let mul = circuit.mul_gate(a_bit, *minus_t_bar);
        let result = circuit.add_gate(b_bit, mul);
        outputs.push(result);
    }

    circuit.output(outputs);
    circuit
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use super::*;
    use crate::{
        bgg::BitToInt,
        poly::{
            dcrt::{params::DCRTPolyParams, poly::DCRTPoly, DCRTPolyUniformSampler},
            enc::rlwe_encrypt,
            sampler::{DistType, PolyUniformSampler},
            Poly, PolyParams,
        },
    };

    #[test]
    fn test_build_composite_circuit_from_public_and_fhe_dec() {
        // 1. Set up parameters
        let params = DCRTPolyParams::default();
        let log_q = params.modulus_bits();
        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sigma = 3.0;

        // 2. Create a simple public circuit that takes 2*log_q inputs and outputs them directly
        let mut public_circuit = PolyCircuit::new();
        {
            let inputs = public_circuit.input(2 * log_q);
            public_circuit.output(inputs[0..2 * log_q].to_vec());
        }

        // 3. Generate a random hardcoded key
        let hardcoded_key = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

        // 4. Generate RLWE ciphertext for the hardcoded key
        let a_rlwe_bar = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);
        let t_bar_matrix = sampler_uniform.sample_uniform(&params, 1, 1, DistType::BitDist);

        let b = rlwe_encrypt(
            &params,
            &sampler_uniform,
            &t_bar_matrix,
            &a_rlwe_bar,
            &hardcoded_key,
            sigma,
        );

        // 5. Build a composite circuit from public circuit and FHE decryption
        let circuit =
            build_composite_circuit_from_public_and_fhe_dec::<DCRTPoly>(public_circuit, log_q);

        // 6. Evaluate the circuit with inputs a_bit_0, b_bit_0, a_bit_1, b_bit_1, ..., -t_bar
        let a_decomposed = a_rlwe_bar.entry(0, 0).decompose_bits(&params);
        let b_decomposed = b.entry(0, 0).decompose_bits(&params);
        let minus_t_bar = -t_bar_matrix.entry(0, 0);
        let mut inputs = Vec::with_capacity(2 * log_q + 1);
        for i in 0..log_q {
            inputs.push(a_decomposed[i].clone());
            inputs.push(b_decomposed[i].clone());
        }
        inputs.push(minus_t_bar.clone());

        assert_eq!(inputs.len(), 2 * log_q + 1);
        let one = DCRTPoly::const_one(&params);
        let outputs = circuit.eval(&params, &one, &inputs);
        assert_eq!(outputs.len(), log_q);

        // 7. Verify the correctness of the output
        for i in 0..log_q {
            let a_bit = a_decomposed[i].clone();
            let b_bit = b_decomposed[i].clone();
            let expected_output = a_bit * minus_t_bar.clone() + b_bit;
            assert_eq!(outputs[i], expected_output);
        }

        // 8. Recompose the output
        let output_ints = DCRTPoly::bits_to_int(&outputs, &params);
        let output_recovered_bits = output_ints.extract_bits_with_threshold(&params);
        assert_eq!(output_recovered_bits, hardcoded_key.entry(0, 0).to_bool_vec());
    }
}
