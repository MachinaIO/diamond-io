use std::collections::HashMap;

use rayon::iter::{IntoParallelIterator, ParallelIterator};
use tracing::info;

use crate::{
    bgg::lut::utils::bgg_encodings_and_input,
    poly::{
        dcrt::{
            sampler::trapdoor::DCRTTrapdoor, DCRTPoly, DCRTPolyMatrix, DCRTPolyParams,
            DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler,
        },
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        PolyMatrix, PolyParams,
    },
};

/// Public Lookup Table
pub struct PublicLut {
    // new public BGG+matrix common for all rows: (n+1)xm, (n+1)xm
    a_lt: (DCRTPolyMatrix, DCRTPolyMatrix),
    // public matrix different for all k: (n+1)xm
    // k => (R_k, L_k)
    lookup_hashmap: HashMap<usize, (DCRTPolyMatrix, DCRTPolyMatrix)>,
    b_l: DCRTPolyMatrix,
}

impl PublicLut {
    pub fn new(
        params: &DCRTPolyParams,
        n: usize,
        f: HashMap<usize, (DCRTPoly, DCRTPoly)>,
        input_size: usize,
        trapdoor: DCRTTrapdoor,
        b_l: DCRTPolyMatrix,
        trap_sampler: DCRTPolyTrapdoorSampler,
    ) -> Self {
        let m = (1 + n) * params.modulus_digits();
        let uni = DCRTPolyUniformSampler::new();
        // public B_L: (n+1)L'xL'(n+1)[logq]
        info!("b_l ({}, {})", b_l.row_size(), b_l.col_size());
        // new public BGG+matrix common for all rows: (n+1)x2m
        let a_lt = uni.sample_uniform(&params, n + 1, 2 * m, DistType::BitDist);
        info!("a_lt ({}, {})", a_lt.row_size(), a_lt.col_size());
        let hashmap_vec: Vec<(usize, (DCRTPolyMatrix, DCRTPolyMatrix))> = (0..f.len())
            .into_par_iter()
            .map(|k| {
                let uni = DCRTPolyUniformSampler::new();
                // public matrix different for all k: (n+1)xm
                let r_k = uni.sample_uniform(params, n + 1, m, DistType::BitDist);
                info!("r_k ({}, {})", r_k.row_size(), r_k.col_size());
                let (x_k, y_k) = f.get(&k).expect("missing f(k)");
                // computing rhs = A_{LT} - (R_k, (- x_k * R_k + y_k * G))
                // todo: check size
                let rhs_rhs =
                    &r_k.concat_columns(&[
                        &(DCRTPolyMatrix::gadget_matrix(params, n + 1) * y_k - &r_k * x_k)
                    ]);
                info!("rhs_rhs ({}, {})", rhs_rhs.row_size(), rhs_rhs.col_size());
                // rhs: (n+1)x2m
                let rhs = &a_lt - rhs_rhs;
                info!("rhs ({}, {})", rhs.row_size(), rhs.col_size());
                // computing target u_1_L' ⊗ rhs: (n+1)L'x2m
                let zeros =
                    DCRTPolyMatrix::zero(params, input_size * rhs.row_size(), rhs.col_size());
                let target = rhs.concat_rows(&[&zeros]);
                info!("target ({}, {})", target.row_size(), target.col_size());
                (k, (r_k, trap_sampler.preimage(params, &trapdoor, &b_l, &target)))
            })
            .collect();
        let a_lt_0 = a_lt.slice(0, n + 1, 0, m);
        let a_lt_1 = a_lt.slice(0, n + 1, m, 2 * m);
        Self { a_lt: (a_lt_0, a_lt_1), lookup_hashmap: hashmap_vec.into_iter().collect(), b_l }
    }

    pub fn evaluate(
        &self,
        params: &DCRTPolyParams,
        n: usize,
        inputs: Vec<usize>,
        p_sigma: f64,
        k: usize,
    ) -> DCRTPolyMatrix {
        let m = (1 + n) * params.modulus_digits();
        let (r_k, l_k) = &self.lookup_hashmap.get(&k).unwrap();
        let (c_x_k, p_x_l) = bgg_encodings_and_input(m, &self.b_l, inputs, n, params, p_sigma);
        // lhs := c_{x_L}G^{-1}(R_k)
        let lhs = c_x_k * r_k.decompose();
        // rhs := p_{x_L}L_k
        let rhs = p_x_l * l_k;
        lhs + rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{poly::Poly, utils::init_tracing};
    use tracing::info;

    #[test]
    fn test_public_lut() {
        init_tracing();
        //  input size 4, dimension 4
        let params = DCRTPolyParams::default();
        // for test purpose T, the number of rows in lookup table is 8,
        let mut f = HashMap::new();
        f.insert(0, (DCRTPoly::const_int(&params, 0), DCRTPoly::const_int(&params, 7)));
        f.insert(1, (DCRTPoly::const_int(&params, 1), DCRTPoly::const_int(&params, 5)));
        // f.insert(2, (DCRTPoly::const_int(&params, 2), DCRTPoly::const_int(&params, 6)));
        // f.insert(3, (DCRTPoly::const_int(&params, 3), DCRTPoly::const_int(&params, 1)));
        // f.insert(4, (DCRTPoly::const_int(&params, 4), DCRTPoly::const_int(&params, 0)));
        // f.insert(5, (DCRTPoly::const_int(&params, 5), DCRTPoly::const_int(&params, 3)));
        // f.insert(6, (DCRTPoly::const_int(&params, 6), DCRTPoly::const_int(&params, 4)));
        // f.insert(7, (DCRTPoly::const_int(&params, 7), DCRTPoly::const_int(&params, 2)));
        let t = f.len();
        let d = 1;
        let inputs = vec![1, 0];
        let input_size = inputs.len();
        let trap_sampler = DCRTPolyTrapdoorSampler::new(&params, 4.578);
        let (trapdoor, b_l) = trap_sampler.trapdoor(&params, (1 + input_size) * (d + 1));
        info!("t:{}, d:{}, input_size:{}", t, d, input_size,);
        let lut = PublicLut::new(&params, d, f, input_size, trapdoor, b_l, trap_sampler);
        let _c_y_k = lut.evaluate(&params, d, inputs, 0.0, 0);
    }
}
