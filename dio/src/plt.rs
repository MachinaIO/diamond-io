use std::collections::HashMap;

use diamond_io::{
    bgg::lut::public_lut::PublicLut,
    poly::{
        Poly,
        dcrt::{
            DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyUniformSampler,
        },
    },
};
use keccak_asm::Keccak256;
use rand::{Rng, rng};

pub fn setup_plt(t_n: usize, params: &DCRTPolyParams, d: usize) -> PublicLut<DCRTPolyMatrix> {
    let mut f = HashMap::new();
    let mut rng = rng();
    for k in 0..t_n {
        let r_val: usize = rng.random_range(0..t_n as usize);
        f.insert(
            k,
            (
                DCRTPoly::from_const_int_lsb(&params, k),
                DCRTPoly::from_const_int_lsb(&params, r_val),
            ),
        );
    }

    let lut = PublicLut::<DCRTPolyMatrix>::new::<
        DCRTPolyUniformSampler,
        DCRTPolyHashSampler<Keccak256>,
    >(params, d, f, rand::random());
    lut
}
