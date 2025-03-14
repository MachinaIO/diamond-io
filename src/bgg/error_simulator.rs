use super::{circuit::PolyCircuit, Evaluable};
use crate::poly::{Poly, PolyElem};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Sub};

impl<P: Poly> PolyCircuit<P> {
    pub fn simulate_error(&self, dim: u32) -> Vec<MPolyCoeffs> {
        let one = ErrorSimulator::default_from_dim(dim);
        let inputs = vec![ErrorSimulator::default_from_dim(dim); self.num_input()];
        let outputs = self.eval(&(), one, &inputs);
        outputs.into_iter().map(|output| output.h_norm).collect_vec()
    }
}

#[derive(Debug, Clone)]
pub struct ErrorSimulator {
    pub h_norm: MPolyCoeffs,
    pub plaintext_norm: BigUint,
    pub dim: u32,
}

impl ErrorSimulator {
    pub fn new(h_norm: MPolyCoeffs, plaintext_norm: BigUint, dim: u32) -> Self {
        Self { h_norm, plaintext_norm, dim }
    }

    pub fn default_from_dim(dim: u32) -> Self {
        let h_norm = MPolyCoeffs::new(vec![BigUint::from(dim)]);
        let plaintext_norm = BigUint::one();
        Self { h_norm, plaintext_norm, dim }
    }
}

impl Add<Self> for ErrorSimulator {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl Add<&Self> for ErrorSimulator {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self {
        Self {
            h_norm: self.h_norm + &rhs.h_norm,
            plaintext_norm: self.plaintext_norm + &rhs.plaintext_norm,
            dim: self.dim,
        }
    }
}

impl Sub<Self> for ErrorSimulator {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl Sub<&Self> for ErrorSimulator {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self {
        Self {
            h_norm: self.h_norm + &rhs.h_norm,
            plaintext_norm: self.plaintext_norm + &rhs.plaintext_norm,
            dim: self.dim,
        }
    }
}

impl Mul<Self> for ErrorSimulator {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self * &rhs
    }
}

impl Mul<&Self> for ErrorSimulator {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self {
        Self {
            h_norm: self.h_norm.right_rotate(self.dim) + rhs.h_norm.clone() * &self.plaintext_norm,
            plaintext_norm: self.plaintext_norm * &rhs.plaintext_norm,
            dim: self.dim,
        }
    }
}

impl<P: Poly> Evaluable<P> for ErrorSimulator {
    type Params = ();
    fn scalar_mul(&self, _: &(), scalar: &P) -> Self {
        let dim = self.dim;
        let scalar_norm =
            scalar.coeffs().iter().fold(BigUint::zero(), |acc, x| acc + x.to_biguint());
        Self {
            h_norm: self.h_norm.right_rotate(dim),
            plaintext_norm: self.plaintext_norm.clone() * scalar_norm,
            dim,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
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

impl Add<&Self> for MPolyCoeffs {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        let self_len = self.0.len();
        let rhs_len = rhs.0.len();
        let max_len = self_len.max(rhs_len);

        let mut result = Vec::with_capacity(max_len);

        for i in 0..max_len {
            let a = if i < self_len { self.0[i].clone() } else { BigUint::zero() };
            let b = if i < rhs_len { rhs.0[i].clone() } else { BigUint::zero() };
            result.push(a + b);
        }

        Self(result)
    }
}

impl Mul<BigUint> for MPolyCoeffs {
    type Output = Self;
    fn mul(self, rhs: BigUint) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&BigUint> for MPolyCoeffs {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::{
        dcrt::{params::DCRTPolyParams, poly::DCRTPoly, sampler::uniform::DCRTPolyUniformSampler},
        sampler::DistType,
    };
    use num_bigint::BigUint;
    use num_traits::Zero;

    fn create_test_error_simulator(
        ring_dim: u32,
        h_norm_values: Vec<u32>,
        plaintext_norm: u32,
    ) -> ErrorSimulator {
        let h_norm = MPolyCoeffs::new(h_norm_values.into_iter().map(BigUint::from).collect());
        let plaintext_norm = BigUint::from(plaintext_norm);
        ErrorSimulator::new(h_norm, plaintext_norm, ring_dim)
    }

    #[test]
    fn test_error_simulator_addition() {
        // Create two ErrorSimulator instances
        let sim1 = create_test_error_simulator(8, vec![10u32], 5);
        let sim2 = create_test_error_simulator(8, vec![20u32], 7);

        // Test addition
        let result = sim1.clone() + sim2.clone();

        // Verify the result
        assert_eq!(result.h_norm.0[0], BigUint::from(30u32)); // 10 + 20
        assert_eq!(result.plaintext_norm, BigUint::from(12u32)); // 5 + 7
        assert_eq!(result.dim, 8);
    }

    #[test]
    fn test_error_simulator_subtraction() {
        // Create two ErrorSimulator instances
        let sim1 = create_test_error_simulator(8, vec![10u32], 5);
        let sim2 = create_test_error_simulator(8, vec![20u32], 7);

        // Test subtraction (which is actually addition in this implementation)
        let result = sim1.clone() - sim2.clone();

        // Verify the result (should be the same as addition)
        assert_eq!(result.h_norm.0[0], BigUint::from(30u32)); // 10 + 20
        assert_eq!(result.plaintext_norm, BigUint::from(12u32)); // 5 + 7
        assert_eq!(result.dim, 8);
    }

    #[test]
    fn test_error_simulator_multiplication() {
        // Create two ErrorSimulator instances
        let sim1 = create_test_error_simulator(8, vec![10u32], 5);
        let sim2 = create_test_error_simulator(8, vec![20u32], 7);

        // Test multiplication
        let result = sim1.clone() * sim2.clone();

        // Verify the result
        // h_norm should be sim1.h_norm.right_rotate(8) + sim2.h_norm * sim1.plaintext_norm

        // Check the length of the h_norm vector
        assert_eq!(result.h_norm.0.len(), 2);

        // First element should be 20 * 5 (from sim2.h_norm * sim1.plaintext_norm)
        assert_eq!(result.h_norm.0[0], BigUint::from(100u32)); // 20 * 5

        // Second element should be 10 * 8 (from right_rotate)
        assert_eq!(result.h_norm.0[1], BigUint::from(80u32));

        assert_eq!(result.plaintext_norm, BigUint::from(35u32)); // 5 * 7
        assert_eq!(result.dim, 8);
    }

    #[test]
    fn test_error_simulator_scalar_multiplication() {
        // Create an ErrorSimulator instance
        let sim = create_test_error_simulator(8, vec![10u32], 5);

        // Create a DCRTPoly for scalar multiplication
        let params = DCRTPolyParams::new(8, 2, 17);
        let sampler = DCRTPolyUniformSampler::new();
        let scalar = sampler.sample_poly(&params, &DistType::FinRingDist);

        // Get the sum of scalar coefficients for verification
        let scalar_norm =
            scalar.coeffs().iter().fold(BigUint::zero(), |acc, x| acc + x.to_biguint());

        // Test scalar multiplication
        let result = sim.scalar_mul(&(), &scalar);

        // Verify the result
        // h_norm should be sim.h_norm.right_rotate(8)
        assert_eq!(result.h_norm.0[0], BigUint::from(0u32));
        assert_eq!(result.h_norm.0[1], BigUint::from(80u32)); // 10 * 8

        // plaintext_norm should be sim.plaintext_norm * sum of scalar coeffs
        assert_eq!(result.plaintext_norm, sim.plaintext_norm * scalar_norm);
        assert_eq!(result.dim, 8);
    }

    #[test]
    fn test_simulate_error() {
        // Create a simple circuit: (input1 + input2) * input3
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let inputs = circuit.input(3);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        let mul_gate = circuit.mul_gate(add_gate, inputs[2]);
        circuit.output(vec![mul_gate]);

        // Simulate error using the circuit
        let error_result = circuit.simulate_error(8);

        // Manually calculate the expected error
        // Create ErrorSimulator instances for inputs
        let input1 = ErrorSimulator::default_from_dim(8);
        let input2 = ErrorSimulator::default_from_dim(8);
        let input3 = ErrorSimulator::default_from_dim(8);

        // Perform the operations manually
        let add_result = input1 + input2;
        let mul_result = add_result * input3;

        // Verify the result
        assert_eq!(error_result.len(), 1);
        assert_eq!(error_result[0], mul_result.h_norm);
    }

    #[test]
    fn test_mpoly_coeffs_operations() {
        // Test MPolyCoeffs addition
        let poly1 = MPolyCoeffs::new(vec![BigUint::from(1u32), BigUint::from(2u32)]);
        let poly2 = MPolyCoeffs::new(vec![BigUint::from(3u32), BigUint::from(4u32)]);

        let result = poly1.clone() + poly2.clone();
        assert_eq!(result.0[0], BigUint::from(4u32)); // 1 + 3
        assert_eq!(result.0[1], BigUint::from(6u32)); // 2 + 4

        // Test MPolyCoeffs right_rotate
        let poly = MPolyCoeffs::new(vec![BigUint::from(5u32), BigUint::from(6u32)]);
        let scale = 3u32;

        let result = poly.right_rotate(scale);
        assert_eq!(result.0[0], BigUint::zero()); // First element is always 0
        assert_eq!(result.0[1], BigUint::from(15u32)); // 5 * 3
        assert_eq!(result.0[2], BigUint::from(18u32)); // 6 * 3

        // Test MPolyCoeffs multiplication by BigUint
        let poly = MPolyCoeffs::new(vec![BigUint::from(7u32), BigUint::from(8u32)]);
        let scalar = BigUint::from(2u32);

        let result = poly * scalar;
        assert_eq!(result.0[0], BigUint::from(14u32)); // 7 * 2
        assert_eq!(result.0[1], BigUint::from(16u32)); // 8 * 2
    }
}
