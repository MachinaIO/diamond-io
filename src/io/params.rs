use crate::{
    bgg::circuit::PolyCircuit,
    poly::{Poly, PolyMatrix, PolyParams},
};

#[derive(Debug, Clone)]
pub struct ObfuscationParams<M: PolyMatrix> {
    pub params: <<M as PolyMatrix>::P as Poly>::Params,
    pub switched_modulus: <<<M as PolyMatrix>::P as Poly>::Params as PolyParams>::Modulus,
    pub input_size: usize,
    pub level_width: usize, // number of bits to be inserted at each level
    pub public_circuit: PolyCircuit,
    /// number of secret key polynomials. This used due to module LWE assumption.
    pub d: usize,
    pub encoding_sigma: f64,
    pub hardcoded_key_sigma: f64,
    pub p_sigma: f64,
    pub trapdoor_sigma: f64,
}
