use super::circuit::{Evaluable, PolyCircuit};
use crate::{
    bgg::lut::public_lut::PublicLut,
    impl_binop_with_refs,
    poly::{dcrt::DCRTPolyMatrix, PolyMatrix},
};
use itertools::Itertools;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul, Sub};

impl<M: PolyMatrix> PolyCircuit<M> {
    pub fn simulate_bgg_norm(
        &self,
        dim: u32,
        base_bits: u32,
        packed_input_norms: Vec<BigUint>,
    ) -> NormBounds
    where
        NormSimulator: Evaluable<Matrix = M>,
    {
        let one = NormSimulator::new(
            MPolyCoeffs::new(vec![BigUint::one()]),
            BigUint::one(),
            dim,
            base_bits,
        );
        let inputs = packed_input_norms
            .into_iter()
            .map(|plaintext_norm| {
                NormSimulator::new(
                    MPolyCoeffs::new(vec![BigUint::one()]),
                    plaintext_norm,
                    dim,
                    base_bits,
                )
            })
            .collect::<Vec<_>>();
        let outputs = self.eval(&(), &one, &inputs, None);
        NormBounds::from_norm_simulators(&outputs)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NormBounds {
    h_norms: Vec<Vec<String>>,
}

impl NormBounds {
    pub fn from_norm_simulators(simulators: &[NormSimulator]) -> Self {
        let h_norms = simulators
            .iter()
            .map(|simulator| {
                let h_norm = simulator.h_norm.0.iter().map(|coeff| coeff.to_string()).collect();
                h_norm
            })
            .collect();
        Self { h_norms }
    }
}

// Note: h_norm and plaintext_norm computed here can be larger than the modulus `q`.
// In such a case, the error after circuit evaluation could be too large.
#[derive(Debug, Clone)]
pub struct NormSimulator {
    pub h_norm: MPolyCoeffs,
    pub plaintext_norm: BigUint,
    pub dim_sqrt: u32,
    pub base: u32,
}

impl NormSimulator {
    pub fn new(h_norm: MPolyCoeffs, plaintext_norm: BigUint, dim: u32, base_bits: u32) -> Self {
        let base = 1 << base_bits;
        let dim_sqrt = (dim as f64).sqrt().ceil() as u32;
        Self { h_norm, plaintext_norm, dim_sqrt, base }
    }
}

impl_binop_with_refs!(NormSimulator => Add::add(self, rhs: &NormSimulator) -> NormSimulator {
    NormSimulator {
        h_norm: &self.h_norm + &rhs.h_norm,
        plaintext_norm: &self.plaintext_norm + &rhs.plaintext_norm,
        dim_sqrt: self.dim_sqrt,
        base: self.base,
    }
});

// Note: norm of the subtraction result is bounded by a sum of the norms of the input matrices,
// i.e., |A-B| < |A| + |B|
impl_binop_with_refs!(NormSimulator => Sub::sub(self, rhs: &NormSimulator) -> NormSimulator {
    NormSimulator {
        h_norm: &self.h_norm + &rhs.h_norm,
        plaintext_norm: &self.plaintext_norm + &rhs.plaintext_norm,
        dim_sqrt: self.dim_sqrt,
        base: self.base,
    }
});

impl_binop_with_refs!(NormSimulator => Mul::mul(self, rhs: &NormSimulator) -> NormSimulator {
    NormSimulator {
        h_norm: self.h_norm.right_rotate(self.dim_sqrt as u64 * (self.base as u64 - 1)) + &rhs.h_norm * &(&self.plaintext_norm * BigUint::from(self.dim_sqrt)),
        plaintext_norm: &self.plaintext_norm * &rhs.plaintext_norm * self.dim_sqrt,
        dim_sqrt: self.dim_sqrt,
        base: self.base,
    }
});

impl Evaluable for NormSimulator {
    type Params = ();
    type Matrix = DCRTPolyMatrix;

    fn rotate(&self, _: &Self::Params, _: usize) -> Self {
        self.clone()
    }

    fn from_digits(_: &Self::Params, one: &Self, digits: &[u32]) -> Self {
        let digit_max = digits.iter().max().unwrap();
        let dim_sqrt = one.dim_sqrt;
        let h_norm = one.h_norm.clone() * BigUint::from(digit_max * dim_sqrt);
        let plaintext_norm = one.plaintext_norm.clone() * BigUint::from(*digit_max);
        Self { h_norm, plaintext_norm, dim_sqrt: one.dim_sqrt, base: one.base }
    }

    fn public_lookup(&self, _plt: &PublicLut<Self::Matrix>, _: Option<Self::Matrix>) -> Self {
        todo!()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MPolyCoeffs(Vec<BigUint>);

impl MPolyCoeffs {
    pub fn new(coeffs: Vec<BigUint>) -> Self {
        Self(coeffs)
    }

    pub fn right_rotate(&self, scale: u64) -> Self {
        let mut coeffs = vec![BigUint::zero()];
        coeffs.extend(self.0.iter().map(|coeff| coeff * scale).collect_vec());
        Self(coeffs)
    }
}

impl_binop_with_refs!(MPolyCoeffs => Add::add(self, rhs: &MPolyCoeffs) -> MPolyCoeffs {
    let self_len = self.0.len();
    let rhs_len = rhs.0.len();
    let max_len = self_len.max(rhs_len);

    let mut result = Vec::with_capacity(max_len);

    for i in 0..max_len {
        let a = if i < self_len { self.0[i].clone() } else { BigUint::zero() };
        let b = if i < rhs_len { rhs.0[i].clone() } else { BigUint::zero() };
        result.push(a + b);
    }

    MPolyCoeffs(result)
});

impl Mul<BigUint> for MPolyCoeffs {
    type Output = Self;
    fn mul(self, rhs: BigUint) -> Self::Output {
        self * &rhs
    }
}

impl Mul<&BigUint> for MPolyCoeffs {
    type Output = Self;
    fn mul(self, rhs: &BigUint) -> Self {
        &self * rhs
    }
}

impl Mul<&BigUint> for &MPolyCoeffs {
    type Output = MPolyCoeffs;
    fn mul(self, rhs: &BigUint) -> MPolyCoeffs {
        MPolyCoeffs(self.0.iter().map(|a| a * rhs).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    fn create_test_error_simulator(
        ring_dim: u32,
        h_norm_values: Vec<u32>,
        plaintext_norm: u32,
    ) -> NormSimulator {
        let h_norm = MPolyCoeffs::new(h_norm_values.into_iter().map(BigUint::from).collect());
        let plaintext_norm = BigUint::from(plaintext_norm);
        NormSimulator::new(h_norm, plaintext_norm, ring_dim, 1)
    }

    #[test]
    fn test_error_simulator_addition() {
        // Create two ErrorSimulator instances
        let sim1 = create_test_error_simulator(16, vec![10u32], 5);
        let sim2 = create_test_error_simulator(16, vec![20u32], 7);

        // Test addition
        let result = sim1 + sim2;

        // Verify the result
        assert_eq!(result.h_norm.0[0], BigUint::from(30u32)); // 10 + 20
        assert_eq!(result.plaintext_norm, BigUint::from(12u32)); // 5 + 7
        assert_eq!(result.dim_sqrt.pow(2), 16);
    }

    #[test]
    fn test_error_simulator_subtraction() {
        // Create two ErrorSimulator instances
        let sim1 = create_test_error_simulator(16, vec![10u32], 5);
        let sim2 = create_test_error_simulator(16, vec![20u32], 7);

        // Test subtraction (which is actually addition in this implementation)
        let result = sim1 - sim2;

        // Verify the result (should be the same as addition)
        assert_eq!(result.h_norm.0[0], BigUint::from(30u32)); // 10 + 20
        assert_eq!(result.plaintext_norm, BigUint::from(12u32)); // 5 + 7
        assert_eq!(result.dim_sqrt.pow(2), 16);
    }

    #[test]
    fn test_error_simulator_multiplication() {
        // Create two ErrorSimulator instances
        let sim1 = create_test_error_simulator(16, vec![10u32], 5);
        let sim2 = create_test_error_simulator(16, vec![20u32], 7);

        // Test multiplication
        let result = sim1 * sim2;

        // Verify the result
        // h_norm should be sim1.h_norm.right_rotate(8) + sim2.h_norm * sim1.plaintext_norm

        // Check the length of the h_norm vector
        assert_eq!(result.h_norm.0.len(), 2);

        // First element should be 20 * 5 * 4 (from sim2.h_norm * sim1.plaintext_norm * dim_sqrt)
        assert_eq!(result.h_norm.0[0], BigUint::from(400u32)); // 20 * 5

        // Second element should be 10 * 4 (from right_rotate)
        assert_eq!(result.h_norm.0[1], BigUint::from(40u32));

        assert_eq!(result.plaintext_norm, BigUint::from(140u32)); // 5 * 7 * 4
        assert_eq!(result.dim_sqrt.pow(2), 16);
    }

    #[test]
    fn test_simulate_bgg_norm() {
        // Create a simple circuit: (input1 + input2) * input3
        let mut circuit: PolyCircuit<DCRTPolyMatrix> = PolyCircuit::new();
        let inputs = circuit.input(3);
        let add_gate = circuit.add_gate(inputs[0], inputs[1]);
        let mul_gate = circuit.mul_gate(add_gate, inputs[2]);
        circuit.output(vec![mul_gate]);

        // Simulate norm using the circuit
        let plaintext_norms = vec![BigUint::from(1u32); 3];
        let norms = circuit.simulate_bgg_norm(16, 1, plaintext_norms);

        // Manually calculate the expected norm
        // Create NormSimulator instances for inputs
        let input1 =
            NormSimulator::new(MPolyCoeffs::new(vec![BigUint::one()]), BigUint::one(), 16, 1);
        let input2 = input1.clone();
        let input3 = input1.clone();

        // Perform the operations manually
        let add_result = input1 + input2;
        let mul_result = add_result * input3;
        let expected = NormBounds::from_norm_simulators(&[mul_result]);
        // Verify the result
        assert_eq!(norms.h_norms.len(), 1);
        assert_eq!(norms, expected);
    }
}
