use crate::{
    bgg::sampler::{BGGEncodingSampler, BGGPublicKeySampler},
    poly::{
        dcrt::{
            sampler::trapdoor::DCRTTrapdoor, DCRTPoly, DCRTPolyHashSampler, DCRTPolyMatrix,
            DCRTPolyParams, DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler,
        },
        sampler::{DistType, PolyHashSampler, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
};
use keccak_asm::Keccak256;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::HashMap;
use tracing::info;

const TAG_R_K: &[u8] = b"TAG_R_K";

/// Public Lookup Table
/// Considering adjusting on diamond-io case.
pub struct PublicLut {
    // public matrix different for all k: (n+1)xm
    // k => (L_k,c_z)
    lookup_hashmap: HashMap<usize, (DCRTPolyMatrix, DCRTPolyMatrix)>,
    r_k_hashkey: [u8; 32],

    /* debugging purpose for correctness */
    #[cfg(feature = "debug")]
    s_x_l: DCRTPolyMatrix,
    #[cfg(feature = "debug")]
    a_lt: DCRTPolyMatrix,
}

impl PublicLut {
    pub fn new(
        params: &DCRTPolyParams,
        n: usize,
        f: HashMap<usize, (DCRTPoly, DCRTPoly)>,
        input_size: usize,
        trapdoor: DCRTTrapdoor,
        b_l: &DCRTPolyMatrix,
        trap_sampler: DCRTPolyTrapdoorSampler,
    ) -> Self {
        let m = (1 + n) * params.modulus_digits();
        let uni = DCRTPolyUniformSampler::new();
        let hash_sampler = DCRTPolyHashSampler::<Keccak256>::new();
        let r_k_hashkey: [u8; 32] = rand::random();
        let t = f.len();
        let r_k_s = hash_sampler.sample_hash(
            params,
            r_k_hashkey,
            TAG_R_K,
            n + 1,
            m * t,
            DistType::FinRingDist,
        );

        // public B_L: (n+1)L'xL'(n+1)[logq]
        info!("b_l ({}, {})", b_l.row_size(), b_l.col_size());
        // new public BGG+matrix common for all rows: (n+1)xm
        let a_lt = uni.sample_uniform(&params, n + 1, m, DistType::BitDist);
        info!("a_lt ({}, {})", a_lt.row_size(), a_lt.col_size());

        /* BGG+ encoding setup */
        let secrets = uni.sample_uniform(&params, 1, n, DistType::BitDist).get_row(0);
        let s_x_l = {
            let minus_one_poly = DCRTPoly::const_minus_one(&params);
            let mut secrets = secrets.to_vec();
            secrets.push(minus_one_poly);
            DCRTPolyMatrix::from_poly_vec_row(&params, secrets)
        };
        info!("s_x_l ({},{})", s_x_l.row_size(), s_x_l.col_size());
        let key: [u8; 32] = rand::random();
        let reveal_plaintexts = vec![false; 2];
        let tag: u64 = rand::random();
        let tag_bytes = tag.to_le_bytes();

        /* Sample c_z(BGG+ encoding), L_k(Preimage) */
        let hashmap_vec: Vec<(usize, (DCRTPolyMatrix, DCRTPolyMatrix))> = (0..t)
            .into_par_iter()
            .map(|k| {
                let (x_k, y_k) = f.get(&k).expect("missing f(k)");
                let uni = DCRTPolyUniformSampler::new();
                let bgg_encoding_sampler = BGGEncodingSampler::new(params, &secrets, uni, 0.0);
                let bgg_pubkey_sampler =
                    BGGPublicKeySampler::<_, DCRTPolyHashSampler<Keccak256>>::new(key, n);
                let pubkeys = bgg_pubkey_sampler.sample(&params, &tag_bytes, &reveal_plaintexts);
                let plaintexts = vec![x_k.clone()];
                let encodings = bgg_encoding_sampler.sample(&params, &pubkeys, &plaintexts);
                assert_eq!(encodings.len(), 2);
                let c_z = encodings[1].clone();
                let a_z = c_z.clone().pubkey.matrix;

                // public matrix different for all k: (n+1)xm
                let r_k = r_k_s.slice_columns(k * m, (k + 1) * m);
                info!("R_k ({}, {})", r_k.row_size(), r_k.col_size());
                info!("A_z ({}, {})", a_z.row_size(), a_z.col_size());

                // computing rhs = A_{LT} - A_{z}G^{-1}(R_k) + x_{k}R_{k} - y_{k}G : (n+1)xm
                let rhs = &a_lt + &(&r_k * x_k) -
                    &(DCRTPolyMatrix::gadget_matrix(params, n + 1) * y_k) -
                    a_z * &r_k.decompose();
                info!("rhs ({}, {})", rhs.row_size(), rhs.col_size());
                // computing target u_1_L' ⊗ rhs: (n+1)L'xm
                let zeros =
                    DCRTPolyMatrix::zero(params, input_size * rhs.row_size(), rhs.col_size());
                let target = rhs.concat_rows(&[&zeros]);
                info!("target ({}, {})", target.row_size(), target.col_size());

                (k, (trap_sampler.preimage(params, &trapdoor, &b_l, &target), c_z.vector))
            })
            .collect();
        #[cfg(feature = "debug")]
        {
            Self { lookup_hashmap: hashmap_vec.into_iter().collect(), r_k_hashkey, s_x_l, a_lt }
        }
        #[cfg(not(feature = "debug"))]
        {
            Self { lookup_hashmap: hashmap_vec.into_iter().collect(), r_k_hashkey }
        }
    }

    pub fn evaluate(
        &self,
        params: &DCRTPolyParams,
        n: usize,
        t: usize,
        p_x_l: DCRTPolyMatrix,
        k: usize,
    ) -> DCRTPolyMatrix {
        let hash_sampler = DCRTPolyHashSampler::<Keccak256>::new();
        let m = (1 + n) * params.modulus_digits();
        let r_k_s = hash_sampler.sample_hash(
            params,
            self.r_k_hashkey,
            TAG_R_K,
            n + 1,
            m * t,
            DistType::FinRingDist,
        );
        let r_k = r_k_s.slice_columns(k * m, (k + 1) * m);
        let (l_k, c_z) = self.lookup_hashmap.get(&k).unwrap();
        let c_lt_k = &p_x_l * l_k;
        let c_y_k = c_z * r_k.decompose() + c_lt_k;

        c_y_k
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{bgg::lut::utils::p_vector_for_inputs, poly::Poly, utils::init_tracing};

    #[test]
    fn test_public_lut() {
        init_tracing();

        /* Setup */
        let params = DCRTPolyParams::default();
        let d = 1;
        let input_size = 2;

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
        let trap_sampler = DCRTPolyTrapdoorSampler::new(&params, 4.578);
        let (trapdoor, b_l) = trap_sampler.trapdoor(&params, (1 + input_size) * (d + 1));
        let t = f.len();

        let lut = PublicLut::new(&params, d, f.clone(), input_size, trapdoor, &b_l, trap_sampler);

        /* Evaluation Step */
        let inputs = vec![1, 0];
        let k = 6;
        assert_eq!(inputs.len(), input_size);
        let p_x_l = p_vector_for_inputs(&b_l, inputs, &params, 0.0, &lut.s_x_l);
        let c_y_k = lut.evaluate(&params, d, t, p_x_l, k);

        /* Information s and y_k is known Only debugging purpose */
        #[cfg(feature = "debug")]
        {
            let (_, y_k) = f.get(&k).expect("missing f(k)");
            let expected_c_y_k =
                &lut.s_x_l * (&lut.a_lt - &(DCRTPolyMatrix::gadget_matrix(&params, d + 1) * y_k));
            debug_assert_eq!(expected_c_y_k, c_y_k);
        }
    }
}
