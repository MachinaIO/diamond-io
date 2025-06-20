//! P to BGG+ conversion

use crate::{
    bgg::{BggEncoding, BggPublicKey},
    poly::{
        dcrt::{DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler},
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        PolyMatrix,
    },
};

pub struct P2BContext {
    /// Single gadget matrix $B_{L+3}$
    pub b_l3: DCRTPolyMatrix,
    pub att_matrix: DCRTPolyMatrix,
    /// One preimage per attribute slot |L|
    pub k_bgg: Vec<DCRTPolyMatrix>,
}

pub fn setup_p2b(params: &DCRTPolyParams, l: usize, d_prime: usize, trap_sigma: f64) -> P2BContext {
    let trap_sampler = DCRTPolyTrapdoorSampler::new(params, trap_sigma);
    let (trap_inv, b_l3) = trap_sampler.trapdoor(params, d_prime);
    let m3 = b_l3.col_size();

    let d_plus_1 = d_prime / 2;
    let uni = DCRTPolyUniformSampler::new();
    let a_prime = uni.sample_uniform(params, d_plus_1, (l + 1) * m3, DistType::BitDist);

    let k_bgg: Vec<DCRTPolyMatrix> = (0..l)
        .map(|i| {
            let a_prime_i = build_a_prime_i(&a_prime, i, m3);
            let gadget_d_plus_1 = DCRTPolyMatrix::gadget_matrix(params, d_plus_1);
            let missing = 2 * m3 - gadget_d_plus_1.col_size();
            let zeros = DCRTPolyMatrix::zero(params, d_plus_1, missing);
            let padded_g = gadget_d_plus_1.concat_columns(&[&zeros]);
            let top = &a_prime_i - &padded_g;
            let zero_rows = DCRTPolyMatrix::zero(params, d_plus_1, 2 * m3);
            let target = top.concat_rows(&[&zero_rows]);

            trap_sampler.preimage(params, &trap_inv, &b_l3, &target)
        })
        .collect();
    P2BContext { b_l3, att_matrix: a_prime, k_bgg }
}

pub(crate) fn build_a_prime_i(full: &DCRTPolyMatrix, i: usize, m: usize) -> DCRTPolyMatrix {
    let start = i * m;
    let stop = (i + 2) * m;
    full.slice(0, full.row_size(), start, stop)
}

pub fn apply_p2b(
    plain_row: &DCRTPolyMatrix,
    ctx: &P2BContext,
    idx: usize,
) -> BggEncoding<DCRTPolyMatrix> {
    let vec = plain_row * &ctx.k_bgg[idx];
    let a_prime_i = build_a_prime_i(&ctx.att_matrix, idx, ctx.b_l3.col_size());
    let pubkey = BggPublicKey { matrix: a_prime_i, reveal_plaintext: false };
    BggEncoding::new(vec, pubkey, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_preimage() {
        let params = DCRTPolyParams::default();
        let l = 1;
        let ctx = setup_p2b(&params, l, 4, 0.0);
        let m3 = ctx.b_l3.col_size();
        let d_plus_1 = 2;
        let a_i = ctx.att_matrix.slice(0, d_plus_1, 0, 2 * m3);
        let gadget_d_plus_1 = DCRTPolyMatrix::gadget_matrix(&params, d_plus_1);
        let pad_cols = 2 * m3 - gadget_d_plus_1.col_size();
        let zero_pad = DCRTPolyMatrix::zero(&params, d_plus_1, pad_cols);
        let padded_g = gadget_d_plus_1.concat_columns(&[&zero_pad]);
        let top = &a_i - &padded_g;
        let target = top.concat_rows(&[&DCRTPolyMatrix::zero(&params, d_plus_1, 2 * m3)]);

        assert_eq!(&ctx.b_l3 * &ctx.k_bgg[0], target);
    }

    #[test]
    fn test_apply_p2b_shape() {
        let params = DCRTPolyParams::default();
        let ctx = setup_p2b(&params, 2, 4, 0.0);
        let m3 = ctx.b_l3.col_size();

        let plain = DCRTPolyMatrix::zero(&params, 1, m3);
        let enc0 = apply_p2b(&plain, &ctx, 0);

        assert_eq!(enc0.vector.row_size(), 1);
        assert_eq!(enc0.vector.col_size(), 2 * m3);

        let d_plus_1 = 2;
        assert_eq!(enc0.pubkey.matrix.row_size(), d_plus_1);
        assert_eq!(enc0.pubkey.matrix.col_size(), 2 * m3);
    }
}
