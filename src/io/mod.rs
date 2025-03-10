pub mod eval;
pub mod obf;
pub mod utils;

use crate::bgg::circuit::PolyCircuit;
use crate::bgg::BggEncoding;
use crate::poly::{matrix::*, polynomial::*, PolyParams};

#[derive(Debug, Clone)]
pub struct Obfuscation<M: PolyMatrix> {
    pub hash_key: [u8; 32],
    pub encodings_init: Vec<BggEncoding<M>>,
    pub p_init: M,
    pub m_preimages: Vec<(M, M)>,
    pub n_preimages: Vec<(M, M)>,
    pub k_preimages: Vec<(M, M)>,
    pub final_preimage: M,
}

#[derive(Debug, Clone)]
pub struct ObfuscationParams<M: PolyMatrix> {
    pub params: <<M as PolyMatrix>::P as Poly>::Params,
    pub switched_modulus: <<<M as PolyMatrix>::P as Poly>::Params as PolyParams>::Modulus,
    pub input_size: usize,
    pub public_circuit: PolyCircuit<M::P>,
    pub error_gauss_sigma: f64,
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use keccak_asm::Keccak256;
    use num_bigint::BigUint;

    use crate::poly::{
        dcrt::{
            DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler,
            DCRTPolyUniformSampler,
        },
        PolyParams,
    };

    use super::{eval::eval_obf, obf::obfuscate, *};

    #[test]
    fn test_io_just_mul_enc_and_bit() {
        let params = DCRTPolyParams::default();
        let log_q = params.modulus_bits();
        let switched_modulus = Arc::new(BigUint::from(3u32));
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
            public_circuit,
            error_gauss_sigma: 0.0,
        };

        let sampler_uniform = DCRTPolyUniformSampler::new();
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let sampler_trapdoor = DCRTPolyTrapdoorSampler::new(2, 0.0);
        let mut rng = rand::rng();
        let obfuscation = obfuscate::<DCRTPolyMatrix, _, _, _, _>(
            obf_params.clone(),
            sampler_uniform,
            sampler_hash,
            sampler_trapdoor,
            &mut rng,
        );
        println!("obfuscated");
        let input = [false];
        let sampler_hash = DCRTPolyHashSampler::<Keccak256>::new([0; 32]);
        let output = eval_obf(obf_params, sampler_hash, obfuscation, &input);
        println!("{:?}", output);
    }
}
