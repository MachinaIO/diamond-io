use std::sync::Arc;

use crate::{
    bgg::circuit::{GateId, PolyCircuit},
    gadgets::big_uint::{BigUintPoly, BigUintPolyContext},
    poly::Poly,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
// ref: https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
// N: the modulus (assumed to be less than 2^64)
// R: 2^{limb_bit_size * num_limbs}
#[derive(Debug, Clone)]
pub struct MontgomeryContext<P: Poly> {
    pub big_uint_ctx: Arc<BigUintPolyContext<P>>,
    pub num_limbs: usize,              // Number of limbs for N
    pub const_n: BigUintPoly<P>,       // N mod R
    pub const_r2: BigUintPoly<P>,      // R^2 mod N
    pub const_n_prime: BigUintPoly<P>, // N' s.t. N' * N = -1 mod R
}

impl<P: Poly> PartialEq for MontgomeryContext<P> {
    fn eq(&self, other: &Self) -> bool {
        self.big_uint_ctx == other.big_uint_ctx &&
            self.const_n.limbs == other.const_n.limbs &&
            self.const_r2.limbs == other.const_r2.limbs &&
            self.const_n_prime.limbs == other.const_n_prime.limbs
    }
}

impl<P: Poly> MontgomeryContext<P> {
    pub fn setup(
        circuit: &mut PolyCircuit<P>,
        params: &P::Params,
        limb_bit_size: usize,
        num_limbs: usize,
        n: u64,
    ) -> Self {
        let big_uint_ctx = Arc::new(BigUintPolyContext::setup(circuit, params, limb_bit_size));

        // Calculate R = 2^(limb_bit_size * num_limbs)
        let r_bits = limb_bit_size * num_limbs;
        let r_big = BigUint::one() << r_bits;
        let n_big = BigUint::from(n);

        // Calculate R^2 mod N
        let r2_big = (&r_big * &r_big) % &n_big;
        let r2 = r2_big.iter_u64_digits().next().unwrap_or(0);

        // Calculate N' such that N * N' ≡ -1 (mod R)
        // This is equivalent to finding the modular inverse of N modulo R, then negating it
        let n_prime_big = Self::calculate_n_prime(&n_big, &r_big);
        let n_prime = n_prime_big.iter_u64_digits().next().unwrap_or(0);
        debug_assert_eq!(
            (&n_big * &n_prime_big) % &r_big,
            &r_big - BigUint::one(),
            "N' calculation failed"
        );

        // Create constant gates
        let const_n = BigUintPoly::const_u64(big_uint_ctx.clone(), circuit, n);
        let const_r2 = BigUintPoly::const_u64(big_uint_ctx.clone(), circuit, r2);
        let const_n_prime = BigUintPoly::const_u64(big_uint_ctx.clone(), circuit, n_prime);

        Self { big_uint_ctx, num_limbs, const_n, const_r2, const_n_prime }
    }

    /// Calculate N' such that N * N' ≡ -1 (mod R)
    /// Using the extended Euclidean algorithm to find the modular inverse
    fn calculate_n_prime(n: &BigUint, r: &BigUint) -> BigUint {
        // We need to find N' such that N * N' ≡ -1 (mod R)
        // This is equivalent to N * N' ≡ R - 1 (mod R)
        // So we find the modular inverse of N modulo R, then multiply by (R-1)

        let n_inv =
            Self::mod_inverse(n, r).expect("N and R must be coprime for Montgomery multiplication");
        let minus_one = r - BigUint::one();
        (n_inv * minus_one) % r
    }

    /// Calculate modular inverse using the extended Euclidean algorithm
    fn mod_inverse(a: &BigUint, m: &BigUint) -> Option<BigUint> {
        if m == &BigUint::one() {
            return Some(BigUint::zero());
        }

        let (mut old_r, mut r) = (a.clone(), m.clone());
        let (mut old_s, mut s) = (BigUint::one(), BigUint::zero());

        while !r.is_zero() {
            let quotient = &old_r / &r;
            let temp_r = &old_r - &quotient * &r;
            old_r = std::mem::replace(&mut r, temp_r);

            let temp_s = if &quotient * &s <= old_s {
                &old_s - &quotient * &s
            } else {
                // Handle underflow by adding m
                m - ((&quotient * &s - &old_s) % m)
            };
            old_s = std::mem::replace(&mut s, temp_s);
        }

        if old_r == BigUint::one() {
            Some(old_s % m)
        } else {
            None // No modular inverse exists
        }
    }
}

#[derive(Debug, Clone)]
pub struct MontgomeryPoly<P: Poly> {
    pub ctx: MontgomeryContext<P>,
    pub value: BigUintPoly<P>,
}

impl<P: Poly> MontgomeryPoly<P> {
    pub fn new(
        circuit: &mut PolyCircuit<P>,
        ctx: MontgomeryContext<P>,
        value: BigUintPoly<P>,
    ) -> Self {
        debug_assert_eq!(value.limbs.len(), ctx.num_limbs, "Value limbs do not match context");
        let r2ed = value.mul(&ctx.const_r2, circuit, None);
        let reduced = montogomery_reduce(&ctx, circuit, &r2ed);
        Self { ctx, value: reduced }
    }

    pub fn mul(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        let muled = self.value.mul(&other.value, circuit, None);
        let reduced = montogomery_reduce(&self.ctx, circuit, &muled);
        Self { ctx: self.ctx.clone(), value: reduced }
    }
}

fn montogomery_reduce<P: Poly>(
    ctx: &MontgomeryContext<P>,
    circuit: &mut PolyCircuit<P>,
    t: &BigUintPoly<P>,
) -> BigUintPoly<P> {
    // u = t * N' mod R
    let u = t.mod_limbs(ctx.num_limbs).mul(
        &ctx.const_n_prime,
        circuit,
        Some(ctx.big_uint_ctx.limb_bit_size * ctx.num_limbs),
    );
    debug_assert_eq!(u.limbs.len(), ctx.num_limbs, "U limbs do not match context");
    // v = (t + m * N) / R
    let v = {
        let muled = u.mul(&ctx.const_n, circuit, None);
        let added = muled.add(&t, circuit);
        added.left_shift(ctx.num_limbs)
    };
    // if v >= N, then return v - N else return v
    let (is_less, diff) = v.less_than(&ctx.const_n, circuit);
    let is_geq = circuit.not_gate(is_less);
    diff.cmux(&v, is_geq, circuit)
}
