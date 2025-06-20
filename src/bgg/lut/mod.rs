pub mod b2p;
pub mod p2b;
pub mod partition;
pub mod utils;

pub use b2p::{apply_b2p, setup_b2p};
pub use p2b::{apply_p2b, setup_p2b};

#[cfg(test)]
mod roundtrip_tests {
    use crate::{
        bgg::lut::{
            apply_b2p, apply_p2b, setup_b2p, setup_p2b, utils::random_bgg_encodings_for_bits,
        },
        poly::{
            dcrt::{DCRTPolyParams, DCRTPolyUniformSampler},
            sampler::{DistType, PolyUniformSampler},
            PolyMatrix,
        },
    };

    #[test]
    fn bgg_to_p_to_bgg() {
        let params = DCRTPolyParams::default();
        let d_prime = 4;
        let l = 1;
        let sigma = 0.0;
        let encs = random_bgg_encodings_for_bits(4, 1, &params);
        assert_eq!(encs.len(), 2);
        let a_one = encs[0].pubkey.matrix.clone();
        let a_attr = encs[1].pubkey.matrix.clone();
        let a_i = a_one.concat_columns(&[&a_attr]);

        /* BGG+ to P */
        let b2p_ctx = setup_b2p(&params, l, d_prime, &[a_i], sigma);
        let uni = DCRTPolyUniformSampler::new();
        let k_rows = b2p_ctx.k_p[0].row_size();
        let p_x_l = uni.sample_uniform(&params, 1, k_rows, DistType::BitDist);
        let plain = apply_b2p((&encs[0], &encs[1]), &b2p_ctx, &p_x_l, 0);

        /* P to BGG+ */
        let p2b_ctx = setup_p2b(&params, l, d_prime, sigma);
        let enc_back = apply_p2b(&plain, &p2b_ctx, 0);
        let original = encs[0].concat_vector(&[encs[1].clone()]);
        let enc_back_prefix = enc_back.vector;

        assert_eq!(enc_back_prefix, original);
    }
}
