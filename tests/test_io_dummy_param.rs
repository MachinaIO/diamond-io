#[cfg(test)]
mod test {
    use diamond_io::{
        bgg::circuit::PolyCircuit,
        io::{obf::obfuscate, params::ObfuscationParams},
        poly::{
            dcrt::{
                DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
                DCRTPolyUniformSampler,
            },
            Poly, PolyParams,
        },
        utils::init_tracing,
    };
    use keccak_asm::Keccak256;
    use num_bigint::BigUint;
    use std::sync::Arc;
    use tracing::info;

    #[test]
    fn test_io_just_mul_enc_and_bit() {
        init_tracing();
        let start_time = std::time::Instant::now();
        let params = DCRTPolyParams::new(4, 2, 17);
        let log_q = params.modulus_bits();
        let switched_modulus = Arc::new(BigUint::from(1u32));
        let mut public_circuit = PolyCircuit::new();
        {
            let inputs = public_circuit.input(log_q + 1);
            let mut outputs = vec![];
            let eval_input = inputs[log_q];
            for enc_input in inputs[0..log_q].iter() {
                let muled = public_circuit.and_gate(*enc_input, eval_input);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        let obf_params = ObfuscationParams {
            params,
            switched_modulus,
            input_size: 1,
            public_circuit: public_circuit.clone(),
            d: 3,
            encoding_sigma: 0.0,
            hardcoded_key_sigma: 0.0,
            p_sigma: 0.0,
        };

        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let sampler_trapdoor = DCRTPolyTrapdoorSampler::new();
        let mut rng = rand::rng();
        let obfuscation = obfuscate::<DCRTPolyMatrix, _, _, _, _>(
            obf_params.clone(),
            sampler_uniform,
            sampler_hash,
            sampler_trapdoor,
            &mut rng,
        );
        let obfuscation_time = start_time.elapsed();
        info!("Time to obfuscate: {:?}", obfuscation_time);

        // let input = [true];
        // let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);

        // #[cfg(feature = "test")]
        // let hardcoded_key = obfuscation
        //     .hardcoded_key
        //     .coeffs()
        //     .iter()
        //     .map(|elem| elem.value() != &BigUint::from(0u8))
        //     .collect::<Vec<_>>();
        // let output = obfuscation.eval(obf_params, sampler_hash, &input);
        // let total_time = start_time.elapsed();
        // info!("{:?}", output);
        // info!("Time for evaluation: {:?}", total_time - obfuscation_time);
        // info!("Total time: {:?}", total_time);
        // #[cfg(feature = "test")]
        // assert_eq!(output, hardcoded_key);
    }
}
