//! BGG+ to P conversion

use crate::{
    bgg::BggEncoding,
    poly::{
        dcrt::{DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler},
        sampler::PolyTrapdoorSampler,
        PolyMatrix,
    },
};

/// Immutable data produced by the obfuscator
pub struct B2PContext {
    /// Stacked $G^{-1}(B_{L+1,*})$
    pub ginv_stacked: DCRTPolyMatrix,
    /// One pre-image per attribute slot (|L| items).
    pub k_p: Vec<DCRTPolyMatrix>,
}

pub fn setup(
    params: &DCRTPolyParams,
    l: usize,
    d_prime: usize,
    a_mats: &[DCRTPolyMatrix],
    trap_sigma: f64,
) -> B2PContext {
    let trap_sampler = DCRTPolyTrapdoorSampler::new(params, trap_sigma);
    let (trap_inv, b_mat) = trap_sampler.trapdoor(params, d_prime);
    let ginv_stacked = stack_ginv_blocks(&b_mat);
    let k_p = (0..l)
        .map(|i| {
            let target = compose_target(&params, &a_mats[i], &ginv_stacked);
            trap_sampler.preimage(&params, &trap_inv, &b_mat, &target)
        })
        .collect();

    B2PContext { ginv_stacked, k_p }
}

pub fn apply(
    pair: (&BggEncoding<DCRTPolyMatrix>, &BggEncoding<DCRTPolyMatrix>),
    ctx: &B2PContext,
    p_x_l: &DCRTPolyMatrix,
    idx: usize,
) -> DCRTPolyMatrix {
    let c_one = &pair.0.vector;
    let c_attr = &pair.1.vector;
    let c_row = c_one.concat_columns(&[&c_attr]);
    let c_p_i = &c_row * &ctx.ginv_stacked;
    let v_p_i = p_x_l * &ctx.k_p[idx];
    &c_p_i - &v_p_i
}

pub(crate) fn stack_ginv_blocks<M>(b_mat: &M) -> M
where
    M: PolyMatrix,
{
    let rows = b_mat.row_size();
    let half = rows / 2;

    let b1 = b_mat.slice_rows(0, half).decompose();
    let b2 = b_mat.slice_rows(half, rows).decompose();

    b1.concat_rows(&[&b2])
}

pub(crate) fn compose_target(
    params: &DCRTPolyParams,
    a_i: &DCRTPolyMatrix,
    ginv_stacked: &DCRTPolyMatrix,
) -> DCRTPolyMatrix {
    assert_eq!(
        a_i.col_size(),
        ginv_stacked.row_size(),
        "A_i must have 2m columns so A_i · ginv_stacked is defined"
    );
    let prod = a_i * ginv_stacked;
    let zeros = DCRTPolyMatrix::zero(params, prod.row_size(), prod.col_size());

    prod.concat_rows(&[&zeros])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bgg::lut::{partition::iter_pairs, utils::random_bgg_encodings_for_bits},
        poly::{
            dcrt::DCRTPolyUniformSampler,
            sampler::{DistType, PolyUniformSampler},
            PolyParams,
        },
    };

    #[test]
    fn test_stack_rows_match_log_base_q() {
        let params = DCRTPolyParams::default();
        let uni = DCRTPolyUniformSampler::new();
        let rows = 4;
        let cols = 3;
        let b_mat = uni.sample_uniform(&params, rows, cols, DistType::BitDist);
        let stacked = stack_ginv_blocks(&b_mat);
        let base_bits = params.base_bits();
        let log_base_q = params.crt_bits().div_ceil(base_bits as usize) * params.crt_depth();

        assert_eq!(stacked.row_size(), rows * log_base_q);
        assert_eq!(stacked.col_size(), cols);
    }

    #[test]
    fn test_compose_target_double_rows() {
        let params = DCRTPolyParams::default();
        let uni = DCRTPolyUniformSampler::new();
        let rows_a = 3;
        let cols_a = 6;
        let a_i = uni.sample_uniform(&params, rows_a, cols_a, DistType::BitDist);
        let ginv = uni.sample_uniform(&params, cols_a, 5, DistType::BitDist);
        let target = compose_target(&params, &a_i, &ginv);

        assert_eq!(target.row_size(), rows_a * 2);
        assert_eq!(target.col_size(), ginv.col_size());
    }

    #[test]
    fn test_setup_outputs_consistent() {
        let params = DCRTPolyParams::default();
        let uni = DCRTPolyUniformSampler::new();
        let d_prime = 4;
        let cols_b = 6;
        let b_mat = uni.sample_uniform(&params, d_prime, cols_b, DistType::BitDist);

        // build G^{-1} stack (rows = 2m)
        let ginv = stack_ginv_blocks(&b_mat);
        let two_m = ginv.row_size();

        //craft A_i with (d+1) × 2m
        let d_plus_1 = d_prime / 2;
        let a_i = uni.sample_uniform(&params, d_plus_1, two_m, DistType::BitDist);
        let target = compose_target(&params, &a_i, &ginv);

        // 2(d+1) × m
        assert_eq!(target.row_size(), d_prime);
        assert_eq!(target.col_size(), cols_b);
    }

    #[test]
    fn test_bgg_to_p_vector() {
        use rayon::iter::ParallelIterator;
        let params = DCRTPolyParams::default();
        let d_prime = 4;
        let trap_sigma = 3.2;
        let ell = params.crt_bits().div_ceil(params.base_bits() as usize) * params.crt_depth();
        let m = d_prime * ell;
        let d_plus_1 = d_prime / 2;
        let uni = DCRTPolyUniformSampler::new();
        let a_i = uni.sample_uniform(&params, d_plus_1, m, DistType::BitDist);

        let ctx = setup(&params, 1, d_prime, &[a_i.clone()], trap_sigma);

        let encs = random_bgg_encodings_for_bits(4, 1, &params);
        let pair = iter_pairs(&encs).find_any(|_| true).expect("at least one attribute");

        let p_cols = ctx.k_p[0].row_size();
        let p_x_l = uni.sample_uniform(&params, 1, p_cols, DistType::BitDist);

        let got = apply(pair, &ctx, &p_x_l, 0);

        let c_row = pair.0.vector.concat_columns(&[&pair.1.vector]);
        let expect = &(&c_row * &ctx.ginv_stacked) - &(p_x_l * &ctx.k_p[0]);

        assert_eq!(got, expect);
    }
}
