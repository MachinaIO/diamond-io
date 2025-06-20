use std::collections::HashMap;

use crate::poly::{
    dcrt::{
        DCRTPoly, DCRTPolyMatrix, DCRTPolyParams, DCRTPolyTrapdoorSampler, DCRTPolyUniformSampler,
    },
    sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
    PolyMatrix,
};

pub struct PublicLut {
    // new public BGG+matrix common for all rows: (n+1)xm, (n+1)xm
    a_lt: (DCRTPolyMatrix, DCRTPolyMatrix),
    // public matrix different for all k: (n+1)xm
    l_k_vec: Vec<DCRTPolyMatrix>,
}

impl PublicLut {
    pub fn new(
        params: &DCRTPolyParams,
        n: usize,
        m: usize,
        t: usize,
        trap_sigma: f64,
        f: HashMap<usize, (DCRTPoly, DCRTPoly)>,
        input_size: usize,
    ) -> Self {
        let uni = DCRTPolyUniformSampler::new();
        // new public BGG+matrix common for all rows: (n+1)x2m
        let a_lt = uni.sample_uniform(&params, n + 1, 2 * m, DistType::BitDist);
        let mut l_k_vec = vec![];
        let trap_sampler = DCRTPolyTrapdoorSampler::new(params, trap_sigma);
        for k in 0..t {
            // public matrix different for all k: (n+1)xm
            let r_k = uni.sample_uniform(&params, n + 1, m, DistType::BitDist);
            // sample L_{k}
            let (trapdoor, b_l) = trap_sampler.trapdoor(&params, n + 1);
            let (x_k, y_k) = f.get(&k).unwrap();
            // computing rhs = A_{LT} - (R_k - x_k * R_k + y_k * G)
            let rhs = &a_lt -
                &(r_k.clone() - (r_k * x_k) +
                    DCRTPolyMatrix::gadget_matrix(&params, n + 1) * y_k);
            // computing u_1_L' ⊗ rhs
            let zeros = DCRTPolyMatrix::zero(&params, input_size * rhs.row_size(), rhs.col_size());
            let target = rhs.concat_rows(&[&zeros]);
            let l_k = trap_sampler.preimage(&params, &trapdoor, &b_l, &target);
            l_k_vec.push(l_k);
        }
        // sample l
        let a_lt_0 = a_lt.slice(0, n + 1, 0, m);
        let a_lt_1 = a_lt.slice(0, n + 1, m, 2 * m);
        Self { a_lt: (a_lt_0, a_lt_1), l_k_vec }
    }
}

#[cfg(test)]
mod tests {
    use crate::poly::{Poly, PolyParams};

    use super::*;

    #[test]
    fn test_public_lut() {
        //  input size 4, dimension 4
        let params = DCRTPolyParams::default();
        // for test purpose T, the number of rows in lookup table is 8,
        let mut f = HashMap::new();
        f.insert(0, (DCRTPoly::const_int(&params, 0), DCRTPoly::const_int(&params, 7)));
        f.insert(1, (DCRTPoly::const_int(&params, 1), DCRTPoly::const_int(&params, 5)));
        f.insert(2, (DCRTPoly::const_int(&params, 2), DCRTPoly::const_int(&params, 6)));
        f.insert(3, (DCRTPoly::const_int(&params, 3), DCRTPoly::const_int(&params, 1)));
        f.insert(4, (DCRTPoly::const_int(&params, 4), DCRTPoly::const_int(&params, 0)));
        f.insert(5, (DCRTPoly::const_int(&params, 5), DCRTPoly::const_int(&params, 3)));
        f.insert(6, (DCRTPoly::const_int(&params, 6), DCRTPoly::const_int(&params, 4)));
        f.insert(7, (DCRTPoly::const_int(&params, 7), DCRTPoly::const_int(&params, 2)));

        let d = 3;
        let m = (d + 1) * params.modulus_digits();
        let lut = PublicLut::new(&params, d, m, 8, 0.0, f, 4);
    }
}
