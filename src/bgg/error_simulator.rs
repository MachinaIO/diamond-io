use super::Evaluable;
use crate::poly::{Poly, PolyElem, PolyParams};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::Zero;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone)]
pub struct ErrorSimulator<P: Poly> {
    pub h_inf: MPolyCoeffs,
    pub plaintext_inf: BigUint,
    pub params: P::Params,
}

impl<P: Poly> ErrorSimulator<P> {
    pub fn new(h_inf: MPolyCoeffs, plaintext_inf: BigUint, params: P::Params) -> Self {
        Self { h_inf, plaintext_inf, params }
    }
}

impl<P: Poly> Add<Self> for ErrorSimulator<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a, P: Poly> Add<&'a Self> for ErrorSimulator<P> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self {
        Self {
            h_inf: self.h_inf + &rhs.h_inf,
            plaintext_inf: self.plaintext_inf + &rhs.plaintext_inf,
            params: self.params,
        }
    }
}

impl<P: Poly> Sub<Self> for ErrorSimulator<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a, P: Poly> Sub<&'a Self> for ErrorSimulator<P> {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self {
        Self {
            h_inf: self.h_inf + &rhs.h_inf,
            plaintext_inf: self.plaintext_inf + &rhs.plaintext_inf,
            params: self.params,
        }
    }
}

impl<P: Poly> Mul<Self> for ErrorSimulator<P> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self * &rhs
    }
}

impl<'a, P: Poly> Mul<&'a Self> for ErrorSimulator<P> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self {
        let dim = self.params.ring_dimension();
        Self {
            h_inf: self.h_inf.right_rotate(dim) + rhs.h_inf.clone() * &self.plaintext_inf,
            plaintext_inf: self.plaintext_inf * &rhs.plaintext_inf,
            params: self.params.clone(),
        }
    }
}

impl<P: Poly> Evaluable<P> for ErrorSimulator<P> {
    fn scalar_mul(&self, _: &<P as Poly>::Params, scalar: &P) -> Self {
        let dim = self.params.ring_dimension();
        let scalar_inf =
            scalar.coeffs().iter().fold(BigUint::zero(), |acc, x| acc + x.to_biguint());
        Self {
            h_inf: self.h_inf.right_rotate(dim),
            plaintext_inf: self.plaintext_inf.clone() * scalar_inf,
            params: self.params.clone(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct MPolyCoeffs(Vec<BigUint>);

impl MPolyCoeffs {
    pub fn new(coeffs: Vec<BigUint>) -> Self {
        Self(coeffs)
    }

    pub fn right_rotate(&self, scale: u32) -> Self {
        let mut coeffs = vec![BigUint::zero()];
        coeffs.extend(self.0.iter().map(|coeff| coeff.clone() * scale).collect_vec());
        Self(coeffs)
    }
}

impl Add<Self> for MPolyCoeffs {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a> Add<&'a Self> for MPolyCoeffs {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        Self(self.0.into_iter().zip(rhs.0.iter()).map(|(a, b)| a + b).collect())
    }
}

impl<'a> Mul<BigUint> for MPolyCoeffs {
    type Output = Self;
    fn mul(self, rhs: BigUint) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a BigUint> for MPolyCoeffs {
    type Output = Self;
    fn mul(self, rhs: &BigUint) -> Self {
        Self(self.0.iter().map(|a| a.clone() * rhs).collect())
    }
}

// #[derive(Debug, Clone)]
// struct MNPoly<E: PolyElem> {
//     coeff_of_degs: BTreeMap<(usize, usize), E>,
// }

// impl<E: PolyElem> MNPoly<E> {
//     pub fn new(monos: Vec<MNMono<E>>) -> Self {
//         Self { monos }
//     }
// }

// impl<'a, E: PolyElem> Add<&'a Self> for MNPoly<E> {
//     type Output = Self;
//     fn add(self, rhs: &Self) -> Self::Output {
//         let mut monos = self.monos.clone();
//         monos.extend(rhs.monos.clone());
//         Self { monos }
//     }
// }

// #[derive(Debug, Clone)]
// struct MNMono<E: PolyElem> {
//     coeff: E,
//     m_degree: usize,
//     n_degree: usize,
// }

// impl<E: PolyElem> MNMono<E> {
//     pub fn new(coeff: E, m_degree: usize, n_degree: usize) -> Self {
//         Self { coeff, m_degree, n_degree }
//     }
// }

// impl<E: PolyElem> Add<Self> for MNMono<E> {
//     type Output = Self;
//     fn add(self, rhs: Self) -> Self::Output {
//         self + &rhs
//     }
// }

// impl<'a, E: PolyElem> Add<&'a Self> for MNMono<E> {
//     type Output = Self;
//     fn add(self, rhs: &Self) -> Self::Output {
//         Self {
//             coeff: self.coeff + &rhs.coeff,
//             m_degree: self.m_degree.max(rhs.m_degree),
//             n_degree: self.n_degree.max(rhs.n_degree),
//         }
//     }
// }

// impl<'a, E: PolyElem> Mul<&'a Self> for MNMono<E> {
//     type Output = Self;
//     fn mul(self, rhs: &Self) -> Self::Output {
//         Self {
//             coeff: self.coeff * &rhs.coeff,
//             m_degree: self.m_degree + rhs.m_degree,
//             n_degree: self.n_degree + rhs.n_degree,
//         }
//     }
// }
