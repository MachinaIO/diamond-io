use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggEncoding,
    },
    poly::{
        dcrt::{
            matrix::base::BaseMatrix, DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix,
            DCRTPolyParams, DCRTPolyUniformSampler,
        },
        sampler::{DistType, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
    utils::create_bit_random_poly,
};
use itertools::Itertools;
use keccak_asm::Keccak256;
use tracing::info;

pub fn random_bgg_encodings_for_bits(
    input_size: usize,
    secret_size: usize,
    params: &DCRTPolyParams,
) -> Vec<BggEncoding<BaseMatrix<DCRTPoly>>> {
    let packed_input_size = input_size.div_ceil(params.ring_dimension().try_into().unwrap());
    // Create samplers
    let key: [u8; 32] = rand::random();
    let bgg_pubkey_sampler =
        BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, secret_size);
    let uniform_sampler = DCRTPolyUniformSampler::new();

    // Generate random tag for sampling
    let tag: u64 = rand::random();
    let tag_bytes = tag.to_le_bytes();
    // Create secret and plaintexts
    let secrets = vec![create_bit_random_poly(&params); secret_size];
    let plaintexts = vec![create_bit_random_poly(&params); packed_input_size];

    // Create random public keys
    let reveal_plaintexts = vec![true; packed_input_size + 1];
    let bgg_encoding_sampler = BGGEncodingSampler::new(params, &secrets, uniform_sampler, 0.0);
    let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
    let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
    encodings
}

pub fn bgg_encodings_and_input(
    m: usize,
    b_l: &DCRTPolyMatrix,
    inputs: Vec<usize>,
    d: usize,
    params: &DCRTPolyParams,
    p_sigma: f64,
) -> DCRTPolyMatrix {
    let input_size = inputs.len();
    let uniform_sampler = DCRTPolyUniformSampler::new();
    // Generate random tag for sampling
    // Create secret and plaintexts
    let s_bars = uniform_sampler.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
    let s_x_l = {
        let minus_one_poly = DCRTPoly::const_minus_one(&params);
        let mut secrets = s_bars.to_vec();
        secrets.push(minus_one_poly);
        DCRTPolyMatrix::from_poly_vec_row(&params, secrets)
    };
    info!("s_x_l ({},{})", s_x_l.row_size(), s_x_l.col_size());
    let t_bar = uniform_sampler.sample_uniform(&params, 1, 1, DistType::BitDist);
    let minus_t_bar = -t_bar.entry(0, 0);
    let one = DCRTPoly::const_one(&params);
    let mut plaintexts = vec![one];
    plaintexts.extend(
        (0..input_size - 1)
            .into_iter()
            .map(|i| DCRTPoly::const_int(params, inputs[i]))
            .collect_vec(),
    );
    plaintexts.push(minus_t_bar);
    // Create random public keys
    let p_x_l_error =
        uniform_sampler.sample_uniform(&params, 1, m, DistType::GaussDist { sigma: p_sigma });
    let encoded_bits = DCRTPolyMatrix::from_poly_vec_row(&params, plaintexts);
    info!("encoded_bits ({},{})", encoded_bits.row_size(), encoded_bits.col_size());
    let s_connect = encoded_bits.tensor(&s_x_l);
    info!("s_connect ({},{})", s_connect.row_size(), s_connect.col_size());
    let s_b = s_connect * b_l;
    let p_x_l = s_b + p_x_l_error;
    p_x_l
}
