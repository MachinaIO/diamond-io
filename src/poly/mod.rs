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
