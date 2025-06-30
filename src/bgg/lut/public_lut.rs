//! Public Lookup

use crate::{
    bgg::sampler::BGGPublicKeySampler,
    poly::{
        sampler::{DistType, PolyHashSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use tracing::info;

const TAG_R_K: &[u8] = b"TAG_R_K";

/// Public Lookup Table
/// Considering adjusting on diamond-io case.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PublicLut<M: PolyMatrix> {
    // public matrix different for all k: (n+1)xm
    // k => (R_k, rhs)
    pub lookup_hashmap: HashMap<usize, (M, M)>,
    r_k_hashkey: [u8; 32],
    // k => (x_k, y_k)
    pub f: HashMap<usize, (M::P, M::P)>,
    pub a_lt: M,
}

impl<M: PolyMatrix> PublicLut<M> {
    pub fn new<SU: PolyUniformSampler<M = M>, SH: PolyHashSampler<[u8; 32], M = M>>(
        params: &<M::P as Poly>::Params,
        d: usize,
        f: HashMap<usize, (M::P, M::P)>,
    ) -> Self {
        let m = (1 + d) * params.modulus_digits();
        let uni = SU::new();
        let hash_sampler = SH::new();
        let r_k_hashkey: [u8; 32] = rand::random();
        let t = f.len();
        let r_k_s = hash_sampler.sample_hash(
            params,
            r_k_hashkey,
            TAG_R_K,
            d + 1,
            m * t,
            DistType::FinRingDist,
        );

        // new public BGG+matrix common for all rows: (n+1)xm
        let a_lt = uni.sample_uniform(params, d + 1, m, DistType::BitDist);
        info!("a_lt ({}, {})", a_lt.row_size(), a_lt.col_size());

        let key: [u8; 32] = rand::random();
        let reveal_plaintexts = vec![false; 2];
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        /* Sample c_z(BGG+ encoding), L_k(Preimage) */
        let hashmap_vec: Vec<(usize, (M, M))> = (0..t)
            .into_par_iter()
            .map(|k| {
                let (x_k, y_k) = f.get(&k).expect("missing f(k)");
                let bgg_pubkey_sampler = BGGPublicKeySampler::<_, SH>::new(key, d);
                let mut pubkeys = bgg_pubkey_sampler.sample(params, &tag_bytes, &reveal_plaintexts);
                assert_eq!(pubkeys.len(), 2);
                let a_z = pubkeys.swap_remove(1).matrix;

                // public matrix different for all k: (n+1)xm
                let r_k = r_k_s.slice_columns(k * m, (k + 1) * m);
                info!("R_k ({}, {})", r_k.row_size(), r_k.col_size());
                info!("A_z ({}, {})", a_z.row_size(), a_z.col_size());

                // computing rhs = A_{LT} - A_{z}G^{-1}(R_k) + x_{k}R_{k} - y_{k}G : (n+1)xm
                let rhs = a_lt.clone() + (r_k.clone() * x_k) -
                    &(M::gadget_matrix(params, d + 1) * y_k) -
                    a_z * r_k.decompose();
                info!("rhs ({}, {})", rhs.row_size(), rhs.col_size());
                // // computing target u_1_L' ⊗ rhs: (n+1)L'xm
                // let zeros = M::zero(params, input_size * rhs.row_size(), rhs.col_size());
                // let target = rhs.concat_rows(&[&zeros]);
                // info!("target ({}, {})", target.row_size(), target.col_size());

                (k, (r_k, rhs))
            })
            .collect();
        Self { lookup_hashmap: hashmap_vec.into_iter().collect(), f, r_k_hashkey, a_lt }
    }
}

#[cfg(test)]
mod tests {
    use keccak_asm::Keccak256;

    use super::*;
    use crate::{
        bgg::lut::utils::p_vector_for_inputs,
        poly::{
            dcrt::{
                DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix, DCRTPolyParams,
                DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler,
            },
            sampler::PolyTrapdoorSampler,
        },
        utils::init_tracing,
    };

    #[test]
    fn test_public_lut() {
        init_tracing();

        /* Setup */
        let params = DCRTPolyParams::default();
        let uni = DCRTPolyUniformSampler::new();
        let d = 1;
        // let input_size = 2;

        /* Lookup mapping k => (x_k, y_k) */
        let mut f = HashMap::new();
        f.insert(0, (DCRTPoly::const_int(&params, 0), DCRTPoly::const_int(&params, 7)));
        f.insert(1, (DCRTPoly::const_int(&params, 1), DCRTPoly::const_int(&params, 5)));
        f.insert(2, (DCRTPoly::const_int(&params, 2), DCRTPoly::const_int(&params, 6)));
        f.insert(3, (DCRTPoly::const_int(&params, 3), DCRTPoly::const_int(&params, 1)));
        f.insert(4, (DCRTPoly::const_int(&params, 4), DCRTPoly::const_int(&params, 0)));
        f.insert(5, (DCRTPoly::const_int(&params, 5), DCRTPoly::const_int(&params, 3)));
        f.insert(6, (DCRTPoly::const_int(&params, 6), DCRTPoly::const_int(&params, 4)));
        f.insert(7, (DCRTPoly::const_int(&params, 7), DCRTPoly::const_int(&params, 2)));

        /* Obfuscation Step */
        // let trap_sampler = DCRTPolyTrapdoorSampler::new(&params, 4.578);
        // let (trapdoor, b_l) = trap_sampler.trapdoor(&params, (1 + input_size) * (d + 1));
        let t = f.len();
        /* BGG+ encoding setup */
        let secrets = uni.sample_uniform(&params, 1, d, DistType::BitDist).get_row(0);
        // in reality there should be input insertion step that updates the secret
        let s_x_l = {
            let minus_one_poly = DCRTPoly::const_minus_one(&params);
            let mut secrets = secrets.to_vec();
            secrets.push(minus_one_poly);
            DCRTPolyMatrix::from_poly_vec_row(&params, secrets)
        };
        info!("s_x_l ({},{})", s_x_l.row_size(), s_x_l.col_size());

        let lut = PublicLut::<DCRTPolyMatrix>::new::<
            DCRTPolyUniformSampler,
            DCRTPolyHashSampler<Keccak256>,
        >(&params, d, f);

        // /* Evaluation Step */
        // let inputs = vec![1, 0];
        // let k = 6;
        // assert_eq!(inputs.len(), input_size);

        // let p_x_l = p_vector_for_inputs(&b_l, inputs, &params, 0.0, &lut.s_x_l);
        // let c_y_k = lut.evaluate(&params, t, p_x_l, k);
        // let (_, y_k) = f.get(&k).expect("missing f(k)");
        // let expected_c_y_k =
        //     &lut.s_x_l * (&lut.a_lt - &(DCRTPolyMatrix::gadget_matrix(&params, d + 1) * y_k));
        // debug_assert_eq!(expected_c_y_k, c_y_k);
    }
}
