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
    use rand::Rng;
    use std::sync::Arc;
    use tracing::info;

    const SIGMA: f64 = 4.578;

    #[test]
    fn test_io_just_mul_enc_and_and_bit_large_width() {
        init_tracing();
        let start_time = std::time::Instant::now();
        let params = DCRTPolyParams::new(4, 2, 17, 10);
        let log_base_q = params.modulus_digits();
        let switched_modulus = Arc::new(BigUint::from(1u32));
        let mut public_circuit = PolyCircuit::new();

        // inputs: BaseDecompose(ct), eval_input
        // outputs: BaseDecompose(ct) AND eval_input, BaseDecompose(ct) AND eval_input
        {
            let inputs = public_circuit.input((2 * log_base_q) + 1);
            let mut outputs = vec![];
            let eval_input = inputs[2 * log_base_q];
            for ct_input in inputs[0..2 * log_base_q].iter() {
                let muled = public_circuit.and_gate(*ct_input, eval_input);
                outputs.push(muled);
            }
            for ct_input in inputs[0..2 * log_base_q].iter() {
                let muled = public_circuit.and_gate(*ct_input, eval_input);
                outputs.push(muled);
            }
            public_circuit.output(outputs);
        }

        let obf_params = ObfuscationParams {
            params: params.clone(),
            switched_modulus,
            input_size: 4,
            level_width_exp: 2,
            public_circuit: public_circuit.clone(),
            d: 3,
            encoding_sigma: 0.0,
            hardcoded_key_sigma: 0.0,
            p_sigma: 0.0,
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
        let input = [bool_in, false, false, false];
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let output =
            obfuscation.eval::<_, DCRTPolyTrapdoorSampler>(obf_params, sampler_hash, &input);
        let total_time = start_time.elapsed();
        info!("Time for evaluation: {:?}", total_time - obfuscation_time);
        info!("Total time: {:?}", total_time);
        let n = output.len() / 2;
        let output_1st_gate = output[..n].to_vec();
        let output_2nd_gate = output[n..].to_vec();
        let input_poly = DCRTPoly::from_const(
            &params,
            &FinRingElem::constant(&params.modulus(), bool_in as u64),
        );
        assert_eq!(output_1st_gate, (hardcoded_key.clone() * input_poly.clone()).to_bool_vec());
        assert_eq!(output_2nd_gate, (hardcoded_key * input_poly).to_bool_vec());
    }
}
