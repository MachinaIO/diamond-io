#![allow(clippy::needless_range_loop)]
#![allow(clippy::suspicious_arithmetic_impl)]

pub mod dcrt;
pub mod element;
pub mod matrix;
pub mod params;
pub mod polynomial;
pub mod sampler;

pub use element::PolyElem;
pub use matrix::PolyMatrix;
pub use params::PolyParams;
pub use polynomial::Poly;

use sampler::{DistType, PolyUniformSampler};
use std::sync::Arc;

pub fn rlwe_encrypt<M, SU>(
    params: &Arc<<<M as PolyMatrix>::P as Poly>::Params>,
    sampler_uniform: &Arc<SU>,
    t: &M,
    a: &M,
    m: &M,
    sigma: f64,
) -> M
where
    M: PolyMatrix,
    SU: PolyUniformSampler<M = M>,
{
    assert!(m.row_size() == 1);
    assert!(m.col_size() == 1);
    assert!(t.row_size() == 1);
    assert!(t.col_size() == 1);
    assert!(a.row_size() == 1);
    assert!(a.col_size() == 1);

    // Sample error from Gaussian distribution
    let e = sampler_uniform.sample_uniform(params.as_ref(), 1, 1, DistType::GaussDist { sigma });

    // Use provided scale or calculate half of q
    let scale = M::P::from_const(params, &<M::P as Poly>::Elem::half_q(&params.modulus()));

    // Compute RLWE encryption: t * a + e - (m * scale)
    t.clone() * a + &e - &(m.clone() * &scale)
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use std::sync::Arc;

    use crate::poly::{
        dcrt::{DCRTPolyMatrix, DCRTPolyParams, DCRTPolyUniformSampler},
        rlwe_encrypt,
        sampler::DistType,
        Poly, PolyMatrix,
    };

    #[test]
    fn test_rlwe_encrypt_decrypt() {
        let params = Arc::new(DCRTPolyParams::default());
        let sampler = Arc::new(DCRTPolyUniformSampler::new());
        let sigma = 3.0;

        // Generate random message bits
        let m = sampler.sample_poly(&params, &DistType::BitDist);

        // Encrypt the message
        let a = sampler.sample_poly(&params, &DistType::BitDist);
        let t = sampler.sample_poly(&params, &DistType::BitDist);

        let m_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![m.clone()]);
        let a_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![a.clone()]);
        let t_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![t.clone()]);
        let b = rlwe_encrypt(&params, &sampler, &t_mat, &a_mat, &m_mat, sigma);

        // Decrypt the ciphertext and recover the message bits
        let recovered = (b - (a_mat * t_mat)).entry(0, 0);
        let recovered_bits = recovered.extract_bits_with_threshold(&params);

        // Verify correctness
        assert_eq!(recovered_bits, m.to_bool_vec());
    }
}
