// pub mod eval;
pub mod obf;
// pub mod serde;
pub mod utils;
use crate::{
    bgg::{circuit::PolyCircuit, BggEncoding},
    poly::{Poly, PolyMatrix, PolyParams},
};

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix> {
    pub hash_key: [u8; 32],
    pub enc_hardcoded_key: M,
    pub encodings_init: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages: Vec<Vec<M>>,
    pub n_preimages: Vec<Vec<M>>,
    pub k_preimages: Vec<Vec<M>>,
    pub final_preimage: M,
    #[cfg(test)]
    pub s_init: M,
    #[cfg(test)]
    pub t_bar: <M as PolyMatrix>::P,
    // #[cfg(test)]
    // pub bs: Vec<Vec<M>>,
    #[cfg(test)]
    pub hardcoded_key: <M as PolyMatrix>::P,
    #[cfg(test)]
    pub final_preimage_target: M,
}

#[derive(Debug, Clone)]
pub struct ObfuscationParams<M: PolyMatrix> {
    pub params: <<M as PolyMatrix>::P as Poly>::Params,
    pub switched_modulus: <<<M as PolyMatrix>::P as Poly>::Params as PolyParams>::Modulus,
    pub input_size: usize,
    pub public_circuit: PolyCircuit,
    pub d: usize,
    pub encoding_sigma: f64,
    pub hardcoded_key_sigma: f64,
    pub p_sigma: f64,
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        bgg::circuit::PolyCircuit,
        io::{obf::obfuscate, ObfuscationParams},
        poly::{
            dcrt::{
                DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
                DCRTPolyUniformSampler,
            },
            PolyParams,
        },
        utils::init_tracing,
    };
    use keccak_asm::Keccak256;
    use num_bigint::BigUint;
    use num_traits::Num;
    use std::sync::Arc;

    #[test]
    fn test_io_just_mul_enc_and_bit() {
        init_tracing();
        let start_time = std::time::Instant::now();
        let params = DCRTPolyParams::default();
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
            params: params.clone(),
            switched_modulus,
            input_size: 1,
            public_circuit: public_circuit.clone(),
            d: 1,
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
        println!("Time to obfuscate: {:?}", obfuscation_time);

        // let input = [true];
        // let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        // let hardcoded_key = obfuscation
        //     .hardcoded_key
        //     .coeffs()
        //     .iter()
        //     .map(|elem| elem.value() != &BigUint::from(0u8))
        //     .collect::<Vec<_>>();
        // let output = obfuscation.eval(obf_params, sampler_hash, &input);
        // let total_time = start_time.elapsed();
        // println!("{:?}", output);
        // println!("Time for evaluation: {:?}", total_time - obfuscation_time);
        // println!("Total time: {:?}", total_time);
        // assert_eq!(output, hardcoded_key);
    }

    // #[test]
    // #[ignore]
    // fn test_io_just_mul_enc_and_bit_real_params() {
    //     init_tracing();
    //     let start_time = std::time::Instant::now();
    //     let params = DCRTPolyParams::new(8192, 7, 51);
    //     let log_q = params.modulus_bits();
    //     let switched_modulus = Arc::new(BigUint::from_str_radix("71671831749689734737838152978190216899892655911508785116799651230841339877765150252187381844012976550000", 10).unwrap());
    //     let mut public_circuit = PolyCircuit::new();
    //     {
    //         let inputs = public_circuit.input(log_q + 1);
    //         let mut outputs = vec![];
    //         let eval_input = inputs[log_q];
    //         for enc_input in inputs[0..log_q].iter() {
    //             let muled = public_circuit.and_gate(*enc_input, eval_input);
    //             outputs.push(muled);
    //         }
    //         public_circuit.output(outputs);
    //     }

    //     let obf_params = ObfuscationParams {
    //         params: params.clone(),
    //         switched_modulus,
    //         input_size: 1,
    //         public_circuit: public_circuit.clone(),
    //         d: 2,
    //         encoding_sigma: 3.477984925390326e-48,
    //         hardcoded_key_sigma: 7.754_896_427_200_485e16,
    //         p_sigma: 1.550677652781115e-169,
    //     };

    //     let sampler_uniform = DCRTPolyUniformSampler::new();
    //     let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
    //     let sampler_trapdoor = DCRTPolyTrapdoorSampler::new();
    //     let mut rng = rand::rng();
    //     let obfuscation = obfuscate::<DCRTPolyMatrix, _, _, _, _>(
    //         obf_params.clone(),
    //         sampler_uniform,
    //         sampler_hash,
    //         sampler_trapdoor,
    //         &mut rng,
    //     );
    //     let obfuscation_time = start_time.elapsed();
    //     println!("Time to obfuscate: {:?}", obfuscation_time);

    //     let input = [true];
    //     let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
    //     let hardcoded_key = obfuscation
    //         .hardcoded_key
    //         .coeffs()
    //         .iter()
    //         .map(|elem| elem.value() != &BigUint::from(0u8))
    //         .collect::<Vec<_>>();
    //     let output = obfuscation.eval(obf_params, sampler_hash, &input);
    //     let total_time = start_time.elapsed();
    //     println!("{:?}", output);
    //     println!("Time for evaluation: {:?}", total_time - obfuscation_time);
    //     println!("Total time: {:?}", total_time);
    //     assert_eq!(output, hardcoded_key);
    // }
}
