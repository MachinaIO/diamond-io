#[cfg(test)]
mod test {
    use diamond_io::{
        bgg::circuit::PolyCircuit,
        io::{obf::obfuscate, params::ObfuscationParams},
        poly::{
            dcrt::{
                DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams,
                DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler, FinRingElem,
            },
            sampler::DistType,
            Poly, PolyElem, PolyParams,
        },
        utils::init_tracing,
    };
    use keccak_asm::Keccak256;
    use num_bigint::BigUint;
    use num_traits::Num;
    use rand::Rng;
    use std::sync::Arc;
    use tracing::info;

    const SIGMA: f64 = 4.578;

    #[test]
    #[ignore]
    fn test_io_just_mul_enc_and_bit_real_params() {
        init_tracing();
        let start_time = std::time::Instant::now();
        let params = DCRTPolyParams::new(8192, 12, 51, 16);
        let log_q = params.modulus_bits();
        let switched_modulus = Arc::new(BigUint::from_str_radix("15829145694278690179872161345257420162248406342477557754500035589945422364945084658438108187108114830000000000000000000000000000000000000000000000000000000000000000000000000000", 10).unwrap());
        let mut public_circuit = PolyCircuit::new();

        // inputs: BITS(ct), eval_input
        // outputs: BITS(ct) AND eval_input
        {
            let inputs = public_circuit.input((2 * log_q) + 1);
            let mut outputs = vec![];
            let eval_input = inputs[2 * log_q];
            for ct_input in inputs[0..2 * log_q].iter() {
                let muled = public_circuit.and_gate(*ct_input, eval_input);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        let obf_params = ObfuscationParams {
            params: params.clone(),
            switched_modulus,
            input_size: 1,
            public_circuit: public_circuit.clone(),
            d: 1,
            encoding_sigma: 7.26009020373352e+51,
            hardcoded_key_sigma: 4.652550537829127e+102,
            p_sigma: 2.9297574882928115,
        };

        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let sampler_trapdoor = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
        let mut rng = rand::rng();
        let hardcoded_key = sampler_uniform.sample_poly(&params, &DistType::BitDist);
        let obfuscation = obfuscate::<DCRTPolyMatrix, _, _, _, _>(
            obf_params.clone(),
            sampler_uniform,
            sampler_hash,
            sampler_trapdoor,
            hardcoded_key.clone(),
            &mut rng,
        );
        let obfuscation_time = start_time.elapsed();
        info!("Time to obfuscate: {:?}", obfuscation_time);

        let bool_in = rng.random::<bool>();
        let input = [bool_in];
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let output =
            obfuscation.eval::<_, DCRTPolyTrapdoorSampler>(obf_params, sampler_hash, &input);
        let total_time = start_time.elapsed();
        info!("hardcoded key {:?}", hardcoded_key.coeffs());
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
        assert_eq!(output, (hardcoded_key * input_poly).to_bool_vec());
    }
}
