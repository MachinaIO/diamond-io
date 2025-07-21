use crate::{
    bgg::{
        sampler::{BGGEncodingSampler, BGGPublicKeySampler},
        BggEncoding,
    },
    poly::{
        dcrt::{
            matrix::base::BaseMatrix, DCRTPoly, DCRTPolyHashSampler, DCRTPolyParams,
            DCRTPolyUniformSampler,
        },
        sampler::PolyUniformSampler,
    },
    utils::{create_bit_random_poly, create_random_poly},
};
use keccak_asm::Keccak256;

pub fn random_bgg_encodings(
    input_size: usize,
    secret_size: usize,
    params: &DCRTPolyParams,
) -> Vec<BggEncoding<BaseMatrix<DCRTPoly>>> {
    // Create samplers
    let key: [u8; 32] = rand::random();
    let bgg_pubkey_sampler =
        BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, secret_size);
    let uniform_sampler = DCRTPolyUniformSampler::new();

    // Generate random tag for sampling
    let tag: u64 = rand::random();
    let tag_bytes = tag.to_le_bytes();
    // Create secret and plaintexts
    let secrets = vec![create_bit_random_poly(params); secret_size];
    let plaintexts = vec![create_random_poly(params); input_size];

    // Create random public keys
    let reveal_plaintexts = vec![true; input_size + 1];
    let bgg_encoding_sampler = BGGEncodingSampler::new(params, &secrets, uniform_sampler, 0.0);
    let pubkeys = bgg_pubkey_sampler.sample(params, &tag_bytes, &reveal_plaintexts);
    bgg_encoding_sampler.sample(params, &pubkeys, &plaintexts)
}
