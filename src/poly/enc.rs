use super::{sampler::PolyUniformSampler, Poly, PolyMatrix};
use crate::poly::{element::PolyElem, sampler::DistType, PolyParams};

pub fn rlwe_encrypt<M, SU>(
    params: &<<M as PolyMatrix>::P as Poly>::Params,
    sampler_uniform: &SU,
    t: M,
    a: &M,
    m: M,
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
    let e = sampler_uniform.sample_uniform(params, 1, 1, DistType::GaussDist { sigma });

    // Use provided scale or calculate half of q
    let scale = M::P::from_const(params, &<M::P as Poly>::Elem::half_q(&params.modulus()));

    // Compute RLWE encryption: t * a + e - (m * scale)
    t * a + &e - &(m * &scale)
}

#[cfg(test)]
#[cfg(feature = "test")]
mod tests {
    use crate::poly::{
        dcrt::{DCRTPolyMatrix, DCRTPolyParams, DCRTPolyUniformSampler},
        enc::rlwe_encrypt,
        sampler::DistType,
        Poly, PolyMatrix,
    };

    #[test]
    fn test_rlwe_encrypt_decrypt() {
        let params = DCRTPolyParams::default();
        let sampler = DCRTPolyUniformSampler::new();
        let sigma = 3.0;

        // Generate random message bits
        let m = sampler.sample_poly(&params, &DistType::BitDist);

        // Encrypt the message
        let a = sampler.sample_poly(&params, &DistType::BitDist);
        let t = sampler.sample_poly(&params, &DistType::BitDist);

        let m_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![m.clone()]);
        let a_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![a.clone()]);
        let t_mat = DCRTPolyMatrix::from_poly_vec_row(&params, vec![t.clone()]);
        let b = rlwe_encrypt(&params, &sampler, t_mat.clone(), &a_mat, m_mat, sigma);

        // Decrypt the ciphertext and recover the message bits
        let recovered = (b - (a_mat * t_mat)).entry(0, 0);
        let recovered_bits = recovered.extract_bits_with_threshold(&params);

        // Verify correctness
        assert_eq!(recovered_bits, m.to_bool_vec());
    }
}
