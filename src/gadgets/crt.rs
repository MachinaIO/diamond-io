use std::sync::Arc;

use crate::{
    bgg::circuit::PolyCircuit,
    gadgets::montgomery::{MontgomeryContext, MontgomeryPoly},
    poly::{Poly, PolyParams},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CrtContext<P: Poly> {
    pub mont_ctxes: Vec<MontgomeryContext<P>>,
    // pub reconstruct_coeffs: Vec<GateId>,
}

impl<P: Poly> CrtContext<P> {
    pub fn setup(circuit: &mut PolyCircuit<P>, params: &P::Params, limb_bit_size: usize) -> Self {
        let (moduli, crt_bits, _crt_depth) = params.to_crt();
        let num_limbs = crt_bits.div_ceil(limb_bit_size);
        let mont_ctxes = moduli
            .iter()
            .map(|modulus| {
                MontgomeryContext::setup(circuit, params, limb_bit_size, num_limbs, *modulus)
            })
            .collect();

        // // Compute CRT reconstruction constants: \tilde{q_i} * q_i^*
        // // where \tilde{q_i} = q / q_i and q_i^* = (\tilde{q_i})^{-1} mod q_i
        // let mut reconstruct_coeffs = Vec::new();

        // for i in 0..moduli.len() {
        //     let qi = moduli[i];

        //     // Compute \tilde{q_i} = q / q_i (product of all other moduli)
        //     let mut q_tilde_i = 1u64;
        //     for j in 0..moduli.len() {
        //         if i != j {
        //             q_tilde_i = (q_tilde_i * moduli[j]) % qi;
        //         }
        //     }

        //     // Compute q_i^* = (\tilde{q_i})^{-1} mod q_i using extended Euclidean algorithm
        //     let qi_star: u64 = mod_inverse(q_tilde_i, qi);

        //     // Compute the reconstruction coefficient: \tilde{q_i} * q_i^* mod q_i
        //     // Since we're working mod q_i, this simplifies to q_i^* mod q_i = 1
        //     // But we need the full coefficient for reconstruction
        //     let mut full_q_tilde_i = 1u64;
        //     for j in 0..moduli.len() {
        //         if i != j {
        //             full_q_tilde_i = full_q_tilde_i.saturating_mul(moduli[j]);
        //         }
        //     }

        //     let coeff = full_q_tilde_i.saturating_mul(qi_star);
        //     let coeff_gate = circuit.add_constant(params, coeff.into());
        //     reconstruct_coeffs.push(coeff_gate);
        // }

        Self { mont_ctxes }
    }
}

pub struct CrtPoly<P: Poly> {
    pub ctx: Arc<CrtContext<P>>,
    pub slots: Vec<MontgomeryPoly<P>>,
}

impl<P: Poly> CrtPoly<P> {
    pub fn new(ctx: Arc<CrtContext<P>>, slots: Vec<MontgomeryPoly<P>>) -> Self {
        Self { ctx, slots }
    }

    pub fn add(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        let new_slots =
            self.slots.iter().zip(other.slots.iter()).map(|(a, b)| a.add(b, circuit)).collect();
        Self::new(self.ctx.clone(), new_slots)
    }

    pub fn sub(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        let new_slots =
            self.slots.iter().zip(other.slots.iter()).map(|(a, b)| a.sub(b, circuit)).collect();
        Self::new(self.ctx.clone(), new_slots)
    }

    pub fn mul(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        let new_slots =
            self.slots.iter().zip(other.slots.iter()).map(|(a, b)| a.mul(b, circuit)).collect();
        Self::new(self.ctx.clone(), new_slots)
    }
}

// Extended Euclidean algorithm to compute modular inverse
// fn mod_inverse(a: u64, m: u64) -> u64 {
//     if m == 1 {
//         return 0;
//     }

//     let (mut a, mut m) = (a as i64, m as i64);
//     let (m0, mut x0, mut x1) = (m, 0i64, 1i64);

//     while a > 1 {
//         let q = a / m;
//         let t = m;
//         m = a % m;
//         a = t;
//         let t = x0;
//         x0 = x1 - q * x0;
//         x1 = t;
//     }

//     if x1 < 0 {
//         x1 += m0;
//     }

//     x1 as u64
// }
