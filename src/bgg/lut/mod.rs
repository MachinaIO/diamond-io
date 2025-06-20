pub mod b2p;
pub mod p2b;
pub mod partition;
pub mod utils;

pub use b2p::{apply_b2p, setup_b2p};
pub use p2b::{apply_p2b, setup_p2b};

#[cfg(test)]
mod roundtrip_tests {
    use rayon::iter::ParallelIterator;

    use crate::{
        bgg::lut::{
            apply_b2p, apply_p2b, partition::iter_pairs, setup_b2p, setup_p2b,
            utils::random_bgg_encodings_for_bits,
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
        let c = random_bgg_encodings_for_bits(4, 1, &params);
        assert_eq!(c.len(), 2);
        for c_i in iter_pairs(&c).collect::<Vec<_>>().into_iter() {
            let combined_c_i = c_i.c_one.concat_vector(&[c_i.c_attr.clone()]);
            let a_i = c_i.a_matrix();

            /* BGG+ to P */
            let b2p_ctx = setup_b2p(&params, l, d_prime, &[a_i], sigma);
            let uni = DCRTPolyUniformSampler::new();
            let k_rows = b2p_ctx.k_p[0].row_size();
            let p_x_l = uni.sample_uniform(&params, 1, k_rows, DistType::BitDist);
            let p_i = apply_b2p(c_i, &b2p_ctx, &p_x_l, 0);

            /* P to BGG+ */
            let p2b_ctx = setup_p2b(&params, l, d_prime, sigma);
            let (enc_one_back, enc_attr_back) = apply_p2b(&p_i, &p2b_ctx, 0);
            let enc_back_vec = enc_one_back.concat_vector(&[enc_attr_back]);

            let enc_back_prefix = enc_back_vec.slice(0, 1, 0, combined_c_i.col_size());

            assert_eq!(enc_back_prefix, combined_c_i);
        }
    }
}
