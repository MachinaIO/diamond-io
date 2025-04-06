use super::{utils::split_int64_mat_alt_to_elems, DCRTTrapdoor};
use crate::{
    parallel_iter,
    poly::{
        dcrt::{
            sampler::{trapdoor::KARNEY_THRESHOLD, DCRTPolyUniformSampler},
            DCRTPoly, DCRTPolyMatrix, DCRTPolyParams,
        },
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        Poly, PolyMatrix, PolyParams,
    },
    utils::{debug_mem, log_mem},
};
use openfhe::ffi::DCRTGaussSampGqArbBase;
use rayon::iter::ParallelIterator;
use std::ops::Range;

const SIGMA: f64 = 4.578;
const SPECTRAL_CONSTANT: f64 = 1.8;

pub struct DCRTPolyTrapdoorSampler {
    sigma: f64,
    base: u32,
    c: f64,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(params: &DCRTPolyParams, sigma: f64) -> Self {
        let base = 1 << params.base_bits();
        let c = (base as f64 + 1.0) * SIGMA;
        Self { sigma, base, c }
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
        let g = DCRTPolyMatrix::gadget_matrix(params, size);
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
        let k = params.modulus_digits();
        let s = SPECTRAL_CONSTANT *
            (self.base as f64 + 1.0) *
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
        log_mem("preimage parameters computed");
        let p_hat = trapdoor.sample_pert_square_mat(
            s,
            self.c,
            self.sigma,
            dgg_large_params,
            peikert,
            target_cols,
        );
        log_mem("p_hat generated");
        let perturbed_syndrome = target - &(public_matrix * &p_hat);
        debug_mem("perturbed_syndrome generated");
        let mut z_hat_mat = DCRTPolyMatrix::zero(params, d * k, target_cols);
        let f = |row_offsets: Range<usize>, col_offsets: Range<usize>| -> Vec<Vec<DCRTPoly>> {
            let nrow = row_offsets.len();
            let ncol = col_offsets.len();
            let perturbed_syndromes = perturbed_syndrome.block_entries(row_offsets, col_offsets);
            let decomposed_results = parallel_iter!(0..nrow)
                .map(|i| {
                    let row_results: Vec<_> = parallel_iter!(0..ncol)
                        .map(|j| {
                            let decomposed = decompose_dcrt_gadget(
                                &perturbed_syndromes[i][j],
                                self.c,
                                params,
                                self.base,
                                self.sigma,
                            );
                            (i, j, decomposed)
                        })
                        .collect();
                    row_results
                })
                .flatten()
                .collect::<Vec<_>>();

            let mut block_matrix = vec![vec![DCRTPoly::const_zero(params); ncol]; k * nrow];
            for (i, j, decomposed) in decomposed_results {
                debug_assert_eq!(decomposed[0].len(), 1);
                for (decomposed_idx, vec) in decomposed.iter().enumerate() {
                    block_matrix[i * k + decomposed_idx][j] = vec[0].clone();
                }
            }
            block_matrix
        };
        z_hat_mat.replace_entries_with_expand(0..d, 0..target_cols, k, 1, f);
        log_mem("z_hat_mat generated");
        let r_z_hat = &trapdoor.r * &z_hat_mat;
        debug_mem("r_z_hat generated");
        let e_z_hat = &trapdoor.e * &z_hat_mat;
        debug_mem("e_z_hat generated");
        let z_hat_former = (p_hat.slice_rows(0, d) + r_z_hat)
            .concat_rows(&[&(p_hat.slice_rows(d, 2 * d) + e_z_hat)]);
        let z_hat_latter = p_hat.slice_rows(2 * d, d * (k + 2)) + z_hat_mat;
        log_mem("z_hat generated");
        z_hat_former.concat_rows(&[&z_hat_latter])
    }
}

pub(crate) fn decompose_dcrt_gadget(
    syndrome: &DCRTPoly,
    c: f64,
    params: &DCRTPolyParams,
    base: u32,
    sigma: f64,
) -> Vec<Vec<DCRTPoly>> {
    let depth = params.crt_depth();
    let z_hat_bbi = parallel_iter!(0..depth)
        .flat_map(|tower_idx| gauss_samp_gq_arb_base(syndrome, c, params, base, sigma, tower_idx))
        .collect::<Vec<_>>();
    split_int64_mat_alt_to_elems(&z_hat_bbi, params)
}

// A function corresponding to lines 260-266 in trapdoor-dcrtpoly.cpp and the `GaussSampGqArbBase`
// function provided by OpenFHE.
pub(crate) fn gauss_samp_gq_arb_base(
    syndrome: &DCRTPoly,
    c: f64,
    params: &DCRTPolyParams,
    base: u32,
    sigma: f64,
    tower_idx: usize,
) -> Vec<Vec<i64>> {
    let n = params.ring_dimension();
    let depth = params.crt_depth();
    let k_res = params.modulus_digits() / depth;
    let result = DCRTGaussSampGqArbBase(
        syndrome.get_poly(),
        c,
        n,
        depth,
        k_res,
        base as i64,
        sigma,
        tower_idx,
    );
    debug_assert_eq!(result.len(), n as usize * k_res);
    // let mut matrix = I64Matrix::new_empty(&I64MatrixParams, k_res, n as usize);
    parallel_iter!(0..k_res)
        .map(|i| {
            parallel_iter!(0..n as usize).map(|j| result[i * n as usize + j]).collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

#[cfg(test)]
#[cfg(feature = "test")]
mod test {
    use super::*;
    use crate::poly::{
        dcrt::{sampler::DCRTPolyUniformSampler, DCRTPolyMatrix, DCRTPolyParams},
        sampler::{DistType, PolyTrapdoorSampler, PolyUniformSampler},
        PolyMatrix, PolyParams,
    };

    const SIGMA: f64 = 4.578;

    #[test]
    fn test_decompose_dcrt_gadget() {
        let params = DCRTPolyParams::default();
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let decomposed = DCRTPolyMatrix::from_poly_vec(
            &params,
            decompose_dcrt_gadget(&target.entry(0, 0), 3.0 * SIGMA, &params, 2, SIGMA),
        );
        let gadget_vec = DCRTPolyMatrix::gadget_matrix(&params, 1);
        assert_eq!(gadget_vec * decomposed, target);
    }

    #[test]
    fn test_decompose_dcrt_gadget_base_8() {
        let params = DCRTPolyParams::new(4, 2, 17, 3);
        let uniform_sampler = DCRTPolyUniformSampler::new();
        let target = uniform_sampler.sample_uniform(&params, 1, 1, DistType::FinRingDist);
        let decomposed = DCRTPolyMatrix::from_poly_vec(
            &params,
            decompose_dcrt_gadget(&target.entry(0, 0), 3.0 * SIGMA, &params, 8, SIGMA),
        );
        let gadget_vec = DCRTPolyMatrix::gadget_matrix(&params, 1);
        assert_eq!(gadget_vec * decomposed, target);
    }

    #[test]
    fn test_trapdoor_generation() {
        let size: usize = 3;
        let params = DCRTPolyParams::default();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);

        let (trapdoor, public_matrix) = trapdoor_sampler.trapdoor(&params, size);

        let expected_rows = size;
        let expected_cols = (&params.modulus_digits() + 2) * size;

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
            let k = params.modulus_digits();
            let identity = DCRTPolyMatrix::identity(&params, size * k, None);
            let trapdoor_matrix = trapdoor.r.concat_rows(&[&trapdoor.e, &identity]);
            public_matrix * trapdoor_matrix
        };
        let gadget_matrix = DCRTPolyMatrix::gadget_matrix(&params, size);
        assert_eq!(muled, gadget_matrix);
    }

    #[test]
    fn test_preimage_generation_square() {
        let params = DCRTPolyParams::default();
        let size = 3;
        let k = params.modulus_digits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
        let k = params.modulus_digits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
        let k = params.modulus_digits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
        let k = params.modulus_digits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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

    #[test]
    fn test_preimage_generation_base_8() {
        let params = DCRTPolyParams::new(4, 2, 17, 3);
        let size = 4;
        let target_cols = 6;
        let k = params.modulus_digits();
        let trapdoor_sampler = DCRTPolyTrapdoorSampler::new(&params, SIGMA);
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
