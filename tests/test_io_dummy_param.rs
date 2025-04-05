#[cfg(test)]
mod test {
    use diamond_io::{
        bgg::circuit::PolyCircuit,
        io::{obf::obfuscate, params::ObfuscationParams, utils::encrypt_rlwe},
        poly::{
            dcrt::{
                DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams,
                DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler, FinRingElem,
            },
            sampler::DistType,
            Poly, PolyElem, PolyMatrix, PolyParams,
        },
        utils::init_tracing,
    };
    use keccak_asm::Keccak256;
    use num_bigint::BigUint;
    use rand::Rng;
    use std::sync::Arc;
    use tracing::info;

    const SIGMA: f64 = 4.578;

    #[test]
    fn test_io_just_mul_enc_and_bit() {
        init_tracing();
        let start_time = std::time::Instant::now();
        let params = DCRTPolyParams::new(4, 2, 17);
        let log_q = params.modulus_bits();
        let switched_modulus = Arc::new(BigUint::from(1u32));

        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let sampler_trapdoor = DCRTPolyTrapdoorSampler::new(SIGMA);

        // 1. Generate RLWE ciphertext (a, b) for input k
        // b = a * t - k * q/2 + e
        let rlwe_encryption_sigma = 3.0;

        // Generate random plaintext bits
        let k = sampler_uniform.sample_poly(&params, &DistType::BitDist);

        // Encrypt the plaintext
        let (a, b, t) = encrypt_rlwe(&params, &sampler_uniform, rlwe_encryption_sigma, &k);

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
        let mut public_circuit = PolyCircuit::new();
        let x_id = public_circuit.input(1)[0];

        let mut public_circuit_outputs = Vec::new();

        let a_bits_vecs = a_bits.iter().map(|poly| poly.to_bool_vec()).collect::<Vec<_>>();
        let b_bits_vecs = b_bits.iter().map(|poly| poly.to_bool_vec()).collect::<Vec<_>>();

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

        // 3. Obfuscate the circuit
        let obf_params = ObfuscationParams {
            params: params.clone(),
            switched_modulus,
            input_size: 1,
            public_circuit: public_circuit.clone(),
            d: 3,
            encoding_sigma: 0.0,
            p_sigma: 0.0,
        };
        let mut rng = rand::rng();

        let t_mat = DCRTPolyMatrix::from_poly_vec_column(&params, vec![t.clone()]);

        let obfuscation = obfuscate::<DCRTPolyMatrix, _, _, _, _>(
            obf_params.clone(),
            &t_mat,
            sampler_uniform,
            sampler_hash,
            sampler_trapdoor,
            &mut rng,
        );
        let obfuscation_time = start_time.elapsed();
        info!("Time to obfuscate: {:?}", obfuscation_time);

        // 4. Evaluate the obfuscated circuit
        let bool_in = rng.random::<bool>();
        let input = [bool_in];
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let output =
            obfuscation.eval::<_, DCRTPolyTrapdoorSampler>(obf_params, sampler_hash, &input);
        let total_time = start_time.elapsed();
        info!("k {:?}", k.coeffs());
        info!("input {:?}", input);
        info!("output {:?}", output);
        info!("Time for evaluation: {:?}", total_time - obfuscation_time);
        info!("Total time: {:?}", total_time);

        #[cfg(feature = "test")]
        let input_poly = DCRTPoly::from_const(
            &params,
            &FinRingElem::constant(&params.modulus(), bool_in as u64),
        );
        #[cfg(feature = "test")]
        assert_eq!(output, (k * input_poly).to_bool_vec());
    }
}
