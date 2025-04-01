use super::{
    trapdoor::{DCRTTrapdoor, KARNEY_THRESHOLD},
    utils::{gen_dcrt_gadget_vector, split_int64_vec_alt_to_elems, split_int64_vec_to_elems},
};
use crate::{
    io::params,
    parallel_iter,
    poly::{
        dcrt::{
            matrix::{i64_matrix::I64MatrixParams, I64Matrix},
            sampler::DCRTPolyUniformSampler,
            DCRTPoly, DCRTPolyMatrix, DCRTPolyParams,
        },
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        PolyMatrix, PolyParams,
    },
    utils::debug_mem,
};
use openfhe::ffi::DCRTGaussSampGqArbBase;
use rayon::iter::ParallelIterator;
use std::ops::Range;

const SIGMA: f64 = 4.578;
const SPECTRAL_CONSTANT: f64 = 1.8;

pub struct DCRTPolyTrapdoorSampler {
    sigma: f64,
    c: f64,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(sigma: f64) -> Self {
        // base = 2
        let c = 3.0 * SIGMA;
        Self { sigma, c }
    }

    fn preimage_square(
        &self,
        params: &DCRTPolyParams,
        trapdoor: &DCRTTrapdoor,
        public_matrix: &DCRTPolyMatrix,
        target: &DCRTPolyMatrix,
        s: f64,
        dgg_large_params: (f64, f64, &[f64]),
        peikert: bool,
    ) -> DCRTPolyMatrix {
        // (d * (k+2)) times d
        let p_hat =
            trapdoor.sample_pert_square_mat(s, self.c, self.sigma, dgg_large_params, peikert);
        println!("p_hat generated");
        let perturbed_syndrome = target.clone() - public_matrix.clone() * &p_hat;
        println!("perturbed_syndrome generated");
        let k = params.modulus_bits();
        let d = public_matrix.row_size();

        let mut z_hat_vecs = Vec::with_capacity(d);
        for i in 0..d {
            let mut z_hat_row_vec = Vec::with_capacity(d);
            for j in 0..d {
                let z_hat = decompose_dcrt_gadget(
                    &perturbed_syndrome.entry(i, j),
                    self.c,
                    params,
                    self.sigma,
                );
                z_hat_row_vec.push(z_hat);
            }
            z_hat_vecs.push(
                z_hat_row_vec[0].concat_columns(&z_hat_row_vec[1..].iter().collect::<Vec<_>>()),
            );
        }
        let z_hat_mat = z_hat_vecs[0].concat_rows(&z_hat_vecs[1..].iter().collect::<Vec<_>>());

        let r_z_hat = trapdoor.r.clone() * &z_hat_mat;
        let e_z_hat = trapdoor.e.clone() * &z_hat_mat;
        let z_hat_former = (p_hat.slice_rows(0, d) + r_z_hat)
            .concat_rows(&[&(p_hat.slice_rows(d, 2 * d) + e_z_hat)]);
        let z_hat_latter = p_hat.slice_rows(2 * d, d * (k + 2)) + z_hat_mat;
        println!("z_hat generated");
        z_hat_former.concat_rows(&[&z_hat_latter])
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = DCRTTrapdoor;

    fn trapdoor(
        &self,
        params: &<<Self::M as crate::poly::PolyMatrix>::P as crate::poly::Poly>::Params,
        size: usize,
    ) -> (Self::Trapdoor, Self::M) {
        let trapdoor = DCRTTrapdoor::new(params, size, self.sigma);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let a_bar = uniform_sampler.sample_uniform(params, size, size, DistType::FinRingDist);
        let g_vec = gen_dcrt_gadget_vector(params);
        let g = g_vec.concat_diag(&vec![&g_vec; size - 1]);
        let a0 = a_bar.concat_columns(&[&DCRTPolyMatrix::identity(params, size, None)]);
        let a1 = g - (a_bar * &trapdoor.r + &trapdoor.e);
        let a = a0.concat_columns(&[&a1]);
        (trapdoor, a)
    }

    fn preimage(
        &self,
        params: &<<Self::M as PolyMatrix>::P as crate::poly::Poly>::Params,
        trapdoor: &Self::Trapdoor,
        public_matrix: &Self::M,
        target: &Self::M,
    ) -> Self::M {
        let d = public_matrix.row_size();
        let target_cols = target.col_size();
        assert_eq!(
            target.row_size(),
            d,
            "Target matrix should have the same number of rows as the public matrix"
        );

        let n = params.ring_dimension() as usize;
        let k = params.modulus_bits();
        let s = SPECTRAL_CONSTANT *
            3.0 *
            SIGMA *
            SIGMA *
            (((d * n * k) as f64).sqrt() + ((2 * n) as f64).sqrt() + 4.7);
        let dgg_large_std = (s * s - self.c * self.c).sqrt();
        let peikert = dgg_large_std < KARNEY_THRESHOLD;
        let (dgg_large_mean, dgg_large_table) = {
            let acc: f64 = 5e-32;
            let m = (-2.0 * acc.ln()).sqrt();
            let fin = (dgg_large_std * m).ceil() as usize;

            let mut m_vals = Vec::with_capacity(fin);
            let variance = 2.0 * dgg_large_std * dgg_large_std;
            let mut cusum = 0.0f64;
            for i in 1..=fin {
                cusum += (-(i as f64 * i as f64) / variance).exp();
                m_vals.push(cusum);
            }
            let m_a = 1.0 / (2.0 * cusum + 1.0);
            for i in 0..fin {
                m_vals[i] *= m_a;
            }
            (m_a, m_vals)
        };
        let dgg_large_params = (dgg_large_mean, dgg_large_std, &dgg_large_table[..]);
        let num_block = target_cols.div_ceil(d);
        debug_mem(format!("preimage before loop processing out of {}", num_block));
        let preimage_blocks = parallel_iter!(0..num_block)
            .map(|i| {
                let start_col = i * d;
                let end_col = (start_col + d).min(target_cols);
                let mut target_block = target.slice(0, d, start_col, end_col);
                let is_padded = end_col - start_col < d;
                if is_padded {
                    let zeros = DCRTPolyMatrix::zero(params, d, start_col + d - end_col);
                    target_block = target_block.concat_columns(&[&zeros]);
                }
                debug_mem(format!("preimage iter : start_col = {}", start_col));
                let mut preimage = self.preimage_square(
                    params,
                    trapdoor,
                    public_matrix,
                    &target_block,
                    s,
                    dgg_large_params,
                    peikert,
                );
                if is_padded {
                    preimage = preimage.slice(0, preimage.row_size(), 0, end_col - start_col);
                }
                preimage
            })
            .collect::<Vec<_>>();
        debug_mem(format!("preimage after loop processing out of {}", num_block));
        preimage_blocks[0].concat_columns(&preimage_blocks[1..].iter().collect::<Vec<_>>())
    }
}

pub(crate) fn decompose_dcrt_gadget(
    syndrome: &DCRTPoly,
    c: f64,
    params: &DCRTPolyParams,
    sigma: f64,
) -> DCRTPolyMatrix {
    let depth = params.crt_depth();
    let z_hat_bbi_blocks = parallel_iter!(0..depth)
        .map(|tower_idx| gauss_samp_gq_arb_base(&syndrome, c, params, sigma, tower_idx))
        .collect::<Vec<_>>();
    println!("z_hat_bbi_blocks generated");
    let z_hat_bbi =
        z_hat_bbi_blocks[0].concat_rows(&z_hat_bbi_blocks[1..].iter().collect::<Vec<_>>());
    split_int64_vec_alt_to_elems(&z_hat_bbi, params)
}

// A function corresponding to lines 260-266 in trapdoor-dcrtpoly.cpp and the `GaussSampGqArbBase`
// function provided by OpenFHE.
pub(crate) fn gauss_samp_gq_arb_base(
    syndrome: &DCRTPoly,
    c: f64,
    params: &DCRTPolyParams,
    sigma: f64,
    tower_idx: usize,
) -> I64Matrix {
    let n = params.ring_dimension();
    let depth = params.crt_depth();
    let k_res = params.modulus_bits() / depth;
    let result =
        DCRTGaussSampGqArbBase(syndrome.get_poly(), c, n, depth, k_res, 2, sigma, tower_idx);
    debug_assert_eq!(result.len(), n as usize * k_res);
    let mut matrix = I64Matrix::zero(&I64MatrixParams, k_res, n as usize);
    let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<i64>> {
        parallel_iter!(row_offsets)
            .map(|i| {
                parallel_iter!(col_offsets.clone())
                    .map(|j| result[i * n as usize + j])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    matrix.replace_entries(0..k_res, 0..n as usize, f);
    matrix
}

#[cfg(test)]
#[cfg(feature = "test")]
mod test {
    use super::*;
    use crate::poly::{
        dcrt::{
            matrix::MmapMatrixElem,
            sampler::{trapdoor::utils::gen_dcrt_gadget_vector, DCRTPolyUniformSampler},
            DCRTPolyMatrix, DCRTPolyParams,
        },
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        PolyMatrix, PolyParams,
    };

    const SIGMA: f64 = 4.578;

    #[test]
    fn test_decompose_dcrt_gadget() {
        let params = DCRTPolyParams::default();
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let decomposed = decompose_dcrt_gadget(&target.entry(0, 0), 3.0 * SIGMA, &params, SIGMA);
        let gadget_vec = gen_dcrt_gadget_vector(&params);
        for i in 0..decomposed.row_size() {
            let entry = decomposed.entry(i, 0);
            println!("{}-th entry: {:?}", i, entry);
        }
        assert_eq!(gadget_vec * decomposed, target);
    }

    #[test]
    fn test_trapdoor_generation() {
        let size: usize = 3;
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(SIGMA);
        let params = DCRTPolyParams::default();

        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        let expected_rows = size;
        let expected_cols = (&params.modulus_bits() + 2) * size;

        assert_eq!(
            public_matrix.row_size(),
            expected_rows,
            "Public matrix should have the correct number of rows"
        );
        assert_eq!(
            public_matrix.col_size(),
            expected_cols,
            "Public matrix should have the correct number of columns"
        );

        // Verify that all entries in the matrix are valid DCRTPolys
        for i in 0..public_matrix.row_size() {
            for j in 0..public_matrix.col_size() {
                let poly = public_matrix.entry(i, j);
                assert!(!poly.get_poly().is_null(), "Matrix entry should be a valid DCRTPoly");
            }
        }

        let muled = {
            let k = params.modulus_bits();
            let identity = DCRTPolyMatrix::identity(&params, size * k, None);
            let trapdoor_matrix = trapdoor.r.concat_rows(&[&trapdoor.e, &identity]);
            public_matrix * trapdoor_matrix
        };
        let gadget_vec = gen_dcrt_gadget_vector(&params);
        let gadget_matrix = gadget_vec.concat_diag(&vec![&gadget_vec; size - 1]);
        assert_eq!(muled, gadget_matrix);
    }

    #[test]
    fn test_preimage_generation_square() {
        let params = DCRTPolyParams::default();
        let size = 3;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(0.0);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, size, size, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

        let expected_rows = size * (k + 2);
        let expected_cols = size;

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns"
        );

        // public_matrix * preimage should be equal to target
        let product = public_matrix * &preimage;
        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }

    #[test]
    fn test_preimage_generation_non_square_target_lt() {
        let params = DCRTPolyParams::default();
        let size = 4;
        let target_cols = 2;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(SIGMA);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols < size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

        let expected_rows = size * (k + 2);
        let expected_cols = target_cols; // Preimage should be sliced to match target columns

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns (sliced to match target)"
        );

        // public_matrix * preimage should be equal to target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }

    #[test]
    fn test_preimage_generation_non_square_target_gt_multiple() {
        let params = DCRTPolyParams::default();
        let size = 4;
        let multiple = 2;
        let target_cols = size * multiple;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(SIGMA);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size
        // target_cols is a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

        let expected_rows = size * (k + 2);
        let expected_cols = target_cols;

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns (equal to target columns)"
        );

        // public_matrix * preimage should be equal to target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }

    #[test]
    fn test_preimage_generation_non_square_target_gt_non_multiple() {
        let params = DCRTPolyParams::default();
        let size = 4;
        let target_cols = 6;
        let k = params.modulus_bits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(SIGMA);
        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        // Create a non-square target matrix (size x target_cols) such that target_cols > size
        // target_cols is not a multiple of size
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target =
            uniform_sampler.sample_uniform(&params, size, target_cols, DistType::FinRingDist);

        let preimage = trapdoor_sampler.preimage(&params, &trapdoor, &public_matrix, &target);

        let expected_rows = size * (k + 2);
        let expected_cols = target_cols;

        assert_eq!(
            preimage.row_size(),
            expected_rows,
            "Preimage matrix should have the correct number of rows"
        );

        assert_eq!(
            preimage.col_size(),
            expected_cols,
            "Preimage matrix should have the correct number of columns (equal to target columns)"
        );

        // public_matrix * preimage should be equal to target
        let product = public_matrix * &preimage;

        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }
}
