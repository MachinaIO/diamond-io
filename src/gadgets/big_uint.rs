use crate::{
    bgg::{
        circuit::{GateId, PolyCircuit},
        lut::public_lut::PublicLut,
    },
    poly::Poly,
};
use std::{collections::HashMap, marker::PhantomData, sync::Arc};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BigUintPolyContext<P: Poly> {
    pub limb_bit_size: usize,
    pub const_zero: GateId,
    pub const_base: GateId,
    pub add_lut_ids: (usize, usize),
    pub mul_lut_ids: (usize, usize),
    _p: PhantomData<P>,
}

impl<P: Poly> BigUintPolyContext<P> {
    pub fn setup(circuit: &mut PolyCircuit<P>, params: &P::Params, limb_bit_size: usize) -> Self {
        let base = 1 << limb_bit_size;
        // assume that base < 2^32
        debug_assert!(limb_bit_size < 32);
        let const_zero = circuit.const_zero_gate();
        let const_base = circuit.const_digits_poly(&[base as u32]);
        let add_luts = Self::setup_split_lut(params, base, 2 * base);
        let mul_luts = Self::setup_split_lut(params, base, base * base);
        let add_lut_ids = (
            circuit.register_public_lookup(add_luts.0),
            circuit.register_public_lookup(add_luts.1),
        );
        let mul_lut_ids = (
            circuit.register_public_lookup(mul_luts.0),
            circuit.register_public_lookup(mul_luts.1),
        );
        Self { limb_bit_size, const_zero, const_base, add_lut_ids, mul_lut_ids, _p: PhantomData }
    }

    fn setup_split_lut(
        params: &P::Params,
        base: usize,
        nrows: usize,
    ) -> (PublicLut<P>, PublicLut<P>) {
        let mut f = HashMap::with_capacity(nrows);
        let mut g = HashMap::with_capacity(nrows);
        for k in 0..nrows {
            let input = P::from_usize_to_constant(params, k);
            let output_f = P::from_usize_to_constant(params, k % base);
            let output_g = P::from_usize_to_constant(params, k / base);
            f.insert(input.clone(), (k, output_f));
            g.insert(input, (k, output_g));
        }
        (PublicLut::new(f), PublicLut::new(g))
    }
}

#[derive(Debug, Clone)]
pub struct BigUintPoly<P: Poly> {
    pub ctx: Arc<BigUintPolyContext<P>>,
    pub limbs: Vec<GateId>,
    _p: PhantomData<P>,
}

impl<P: Poly> BigUintPoly<P> {
    pub fn new(ctx: Arc<BigUintPolyContext<P>>, limbs: Vec<GateId>) -> Self {
        Self { ctx, limbs, _p: PhantomData }
    }

    pub fn bit_size(&self) -> usize {
        self.limbs.len() * self.ctx.limb_bit_size
    }

    pub fn zero(ctx: Arc<BigUintPolyContext<P>>, bit_size: usize) -> Self {
        debug_assert_eq!(bit_size % ctx.limb_bit_size, 0);
        let limb_len = bit_size / ctx.limb_bit_size;
        let limbs = vec![ctx.const_zero; limb_len];
        Self { ctx, limbs, _p: PhantomData }
    }

    pub fn const_u64(
        ctx: Arc<BigUintPolyContext<P>>,
        circuit: &mut PolyCircuit<P>,
        value: u64,
    ) -> Self {
        let mut limbs = vec![];
        let mut remaining_value = value;
        let base = 1u64 << ctx.limb_bit_size;
        while remaining_value > 0 {
            limbs.push(circuit.const_digits_poly(&[(remaining_value % base) as u32]));
            remaining_value /= base;
        }

        debug_assert_eq!(remaining_value, 0);
        Self { ctx, limbs, _p: PhantomData }
    }

    pub fn extend_size(&self, new_bit_size: usize) -> Self {
        debug_assert!(new_bit_size >= self.bit_size());
        debug_assert_eq!(new_bit_size % self.ctx.limb_bit_size, 0);
        let limb_len = new_bit_size / self.ctx.limb_bit_size + 1;
        let mut limbs = self.limbs.clone();
        limbs.resize(limb_len, self.ctx.const_zero);
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }

    pub fn print(&self, prefix: String, circuit: &mut PolyCircuit<P>) {
        for (idx, limb) in self.limbs.iter().enumerate() {
            circuit.print(*limb, format!("{prefix}, {idx}-th limb"));
        }
    }

    pub fn add(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        let (max_num_limbs, a, b) = if self.limbs.len() >= other.limbs.len() {
            (self.limbs.len(), &self.limbs, &other.limbs)
        } else {
            (other.limbs.len(), &other.limbs, &self.limbs)
        };
        let mut limbs = Vec::with_capacity(max_num_limbs + 1);
        let mut carry = circuit.const_zero_gate();
        for i in 0..max_num_limbs {
            let sum = if i >= b.len() {
                circuit.add_gate(a[i], carry)
            } else {
                let tmp = circuit.add_gate(a[i], b[i]);
                circuit.add_gate(tmp, carry)
            };
            let sum_l = circuit.public_lookup_gate(sum, self.ctx.add_lut_ids.0);
            let sum_h = circuit.public_lookup_gate(sum, self.ctx.add_lut_ids.1);
            limbs.push(sum_l);
            carry = sum_h;
        }
        limbs.push(carry);
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }

    pub fn less_than(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> (GateId, Self) {
        debug_assert_eq!(self.limbs.len(), other.limbs.len());
        debug_assert_eq!(self.ctx, other.ctx);
        let mut limbs = Vec::with_capacity(self.limbs.len());
        let mut borrow = circuit.const_zero_gate();
        let one = circuit.const_one_gate();
        for i in 0..self.limbs.len() {
            let tmp0 = circuit.add_gate(self.limbs[i], self.ctx.const_base);
            let tmp1 = circuit.sub_gate(tmp0, other.limbs[i]);
            let diff = circuit.sub_gate(tmp1, borrow);
            let diff_l = circuit.public_lookup_gate(diff, self.ctx.add_lut_ids.0);
            let diff_h = circuit.public_lookup_gate(diff, self.ctx.add_lut_ids.1);
            limbs.push(diff_l);
            borrow = circuit.sub_gate(one, diff_h);
        }
        (borrow, Self { ctx: self.ctx.clone(), limbs, _p: PhantomData })
    }

    pub fn mul(
        &self,
        other: &Self,
        circuit: &mut PolyCircuit<P>,
        max_bit_size: Option<usize>,
    ) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        let max_bit_size = max_bit_size.unwrap_or(self.bit_size() + other.bit_size());
        debug_assert!(max_bit_size % self.ctx.limb_bit_size == 0);
        let max_limbs = max_bit_size / self.ctx.limb_bit_size;
        let zero = circuit.const_zero_gate();
        let mut add_limbs = vec![vec![]; max_limbs];
        for i in 0..self.limbs.len() {
            for j in 0..other.limbs.len() {
                if i + j >= max_limbs {
                    continue; // skip if next index exceeds max limbs
                }
                let mul = circuit.mul_gate(self.limbs[i], other.limbs[j]);
                let mul_l = circuit.public_lookup_gate(mul, self.ctx.mul_lut_ids.0);
                add_limbs[i + j].push(mul_l);
                if i + j + 1 >= max_limbs {
                    continue; // skip if next index exceeds max limbs
                }
                let mul_h = circuit.public_lookup_gate(mul, self.ctx.mul_lut_ids.1);
                add_limbs[i + j + 1].push(mul_h);
            }
        }
        let mut limbs = vec![zero; max_limbs];
        for i in 0..add_limbs.len().min(max_limbs) {
            let add_limb = &add_limbs[i];
            if add_limb.is_empty() {
                continue; // skip if no additions for this limb
            }
            let mut carry = circuit.const_zero_gate();
            let mut sum_l = add_limb[0];
            for limb in add_limb.iter().skip(1) {
                let sum = circuit.add_gate(sum_l, *limb);
                sum_l = circuit.public_lookup_gate(sum, self.ctx.add_lut_ids.0);
                let sum_h = circuit.public_lookup_gate(sum, self.ctx.add_lut_ids.1);
                carry = circuit.add_gate(carry, sum_h);
            }
            limbs[i] = sum_l;
            if i + 1 < max_limbs {
                add_limbs[i + 1].push(carry);
            }
        }
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }

    pub fn left_shift(&self, shift: usize) -> Self {
        debug_assert!(shift < self.limbs.len());
        let limbs = self.limbs[shift..].to_vec();
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }

    pub fn mod_limbs(&self, num_limbs: usize) -> Self {
        debug_assert!(num_limbs <= self.limbs.len());
        let limbs = self.limbs[0..num_limbs].to_vec();
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }

    // return self if selector is 1, other if selector is 0
    pub fn cmux(&self, other: &Self, selector: GateId, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.ctx, other.ctx);
        debug_assert_eq!(self.limbs.len(), other.limbs.len());
        let mut limbs = Vec::with_capacity(self.limbs.len());
        let not = circuit.not_gate(selector);
        for i in 0..self.limbs.len() {
            let case1 = circuit.mul_gate(self.limbs[i], selector);
            let case2 = circuit.mul_gate(other.limbs[i], not);
            let cmuxed = circuit.add_gate(case1, case2);
            limbs.push(cmuxed);
        }
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        bgg::circuit::{eval::PolyPltEvaluator, PolyCircuit},
        poly::dcrt::{params::DCRTPolyParams, poly::DCRTPoly},
    };
    use std::sync::Arc;

    const INPUT_BIT_SIZE: usize = 20;
    const LIMB_BIT_SIZE: usize = 5;
    const LIMB_LEN: usize = INPUT_BIT_SIZE / LIMB_BIT_SIZE;

    fn create_test_context(
        circuit: &mut PolyCircuit<DCRTPoly>,
        total_limbs: usize,
    ) -> (Vec<GateId>, DCRTPolyParams, Arc<BigUintPolyContext<DCRTPoly>>) {
        let params = DCRTPolyParams::default();
        let inputs = circuit.input(total_limbs);
        let ctx = Arc::new(BigUintPolyContext::setup(circuit, &params, LIMB_BIT_SIZE));
        (inputs, params, ctx)
    }

    fn create_test_biguint_from_value(
        ctx: Arc<BigUintPolyContext<DCRTPoly>>,
        params: &DCRTPolyParams,
        value: u32,
    ) -> Vec<DCRTPoly> {
        let limb_len = LIMB_LEN;
        let mut limbs = Vec::with_capacity(limb_len);
        let mut remaining_value = value;
        let base = 1u32 << ctx.limb_bit_size;

        for _ in 0..limb_len {
            let limb_value = remaining_value % base;
            limbs.push(DCRTPoly::from_usize_to_constant(&params, limb_value as usize));
            remaining_value /= base;
        }

        debug_assert_eq!(remaining_value, 0);
        limbs
    }

    #[test]
    fn test_biguint_add() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        let result = big_a.add(&big_b, &mut circuit);
        circuit.output(result.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 15);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 20);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        let mut expected_sum = 15u32 + 20u32;
        let mut expected_limbs = vec![0; LIMB_LEN + 1];
        for i in 0..LIMB_LEN + 1 {
            if expected_sum == 0 {
                break;
            }
            expected_limbs[i] = expected_sum % (1u32 << ctx.limb_bit_size);
            expected_sum /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), 5);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_add_with_carry() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        let result = big_a.add(&big_b, &mut circuit);
        circuit.output(result.limbs.clone());

        // Use values that will cause carry with 20-bit input size (4 limbs of 5 bits each)
        let a = create_test_biguint_from_value(ctx.clone(), &params, 1048575); // 2^20 - 1
        let b = create_test_biguint_from_value(ctx.clone(), &params, 1);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        let mut expected_sum = 1048575u32 + 1u32; // This will be 1048576 = 2^20
        let mut expected_limbs = vec![0; LIMB_LEN + 1];
        for i in 0..LIMB_LEN + 1 {
            if expected_sum == 0 {
                break;
            }
            expected_limbs[i] = expected_sum % (1u32 << ctx.limb_bit_size);
            expected_sum /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), LIMB_LEN + 1);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_less_than_smaller() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        let (lt_result, diff) = big_a.less_than(&big_b, &mut circuit);
        let mut output_gates = vec![lt_result];
        output_gates.extend(diff.limbs.clone());
        circuit.output(output_gates);

        // a < b (500 < 1000), so less_than should return 1 (true)
        let a = create_test_biguint_from_value(ctx.clone(), &params, 500);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 1000);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        assert_eq!(eval_result.len(), 1 + LIMB_LEN);

        let lt_coeffs = eval_result[0].coeffs();
        assert_eq!(*lt_coeffs[0].value(), 1u32.into());

        let mut expected_diff = 500u32 + (1u32 << INPUT_BIT_SIZE) - 1000u32;
        let mut expected_limbs = vec![0; LIMB_LEN];
        for i in 0..LIMB_LEN {
            if expected_diff == 0 {
                break;
            }
            expected_limbs[i] = expected_diff % (1u32 << ctx.limb_bit_size);
            expected_diff /= 1u32 << ctx.limb_bit_size;
        }

        for i in 1..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i - 1].into());
        }
    }

    #[test]
    fn test_biguint_less_than_equal() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        let (lt_result, diff) = big_a.less_than(&big_b, &mut circuit);
        let mut output_gates = vec![lt_result];
        output_gates.extend(diff.limbs.clone());
        circuit.output(output_gates);

        // a == b (12345 == 12345), so less_than should return 0 (false)
        let a = create_test_biguint_from_value(ctx.clone(), &params, 12345);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 12345);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        assert_eq!(eval_result.len(), 1 + LIMB_LEN);

        let lt_coeffs = eval_result[0].coeffs();
        assert_eq!(*lt_coeffs[0].value(), 0u32.into());

        let mut expected_diff = 1u32 << INPUT_BIT_SIZE; // base
        let mut expected_limbs = vec![0; LIMB_LEN];
        for i in 0..LIMB_LEN {
            if expected_diff == 0 {
                break;
            }
            expected_limbs[i] = expected_diff % (1u32 << ctx.limb_bit_size);
            expected_diff /= 1u32 << ctx.limb_bit_size;
        }

        for i in 1..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i - 1].into());
        }
    }

    #[test]
    fn test_biguint_less_than_greater() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        let (lt_result, diff) = big_a.less_than(&big_b, &mut circuit);
        let mut output_gates = vec![lt_result];
        output_gates.extend(diff.limbs.clone());
        circuit.output(output_gates);

        // a > b (1000 > 500), so less_than should return 0 (false)
        let a = create_test_biguint_from_value(ctx.clone(), &params, 1000);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 500);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        assert_eq!(eval_result.len(), 1 + LIMB_LEN);

        let lt_coeffs = eval_result[0].coeffs();
        assert_eq!(*lt_coeffs[0].value(), 0u32.into());

        let mut expected_diff = 1000u32 - 500u32;
        let mut expected_limbs = vec![0; LIMB_LEN];
        for i in 0..LIMB_LEN {
            if expected_diff == 0 {
                break;
            }
            expected_limbs[i] = expected_diff % (1u32 << ctx.limb_bit_size);
            expected_diff /= 1u32 << ctx.limb_bit_size;
        }

        for i in 1..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i - 1].into());
        }
    }

    #[test]
    fn test_biguint_mul_simple() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        let result = big_a.mul(&big_b, &mut circuit, None);
        circuit.output(result.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 123);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 456);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        let mut expected_product = 123u32 * 456u32;
        let output_limb_len = 40 / LIMB_BIT_SIZE; // 8 limbs for 40-bit output
        let mut expected_limbs = vec![0; output_limb_len];
        for i in 0..output_limb_len {
            if expected_product == 0 {
                break;
            }
            expected_limbs[i] = expected_product % (1u32 << ctx.limb_bit_size);
            expected_product /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), output_limb_len);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_mul_with_overflow() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..].to_vec());
        // Use values that will cause overflow with 20-bit input size
        let result = big_a.mul(&big_b, &mut circuit, Some(40));
        circuit.output(result.limbs.clone());

        // Use larger values that will produce overflow across multiple limbs
        let a = create_test_biguint_from_value(ctx.clone(), &params, 1023); // near max for 10 bits
        let b = create_test_biguint_from_value(ctx.clone(), &params, 1023);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b].concat(),
            Some(plt_evaluator),
        );

        let mut expected_product = 1023u32 * 1023u32;
        let output_limb_len = 40 / LIMB_BIT_SIZE; // 8 limbs for 40-bit output
        let mut expected_limbs = vec![0; output_limb_len];
        for i in 0..output_limb_len {
            if expected_product == 0 {
                break;
            }
            expected_limbs[i] = expected_product % (1u32 << ctx.limb_bit_size);
            expected_product /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), output_limb_len);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_zero() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (_, params, ctx) = create_test_context(&mut circuit, LIMB_LEN);

        let zero = BigUintPoly::zero(ctx.clone(), INPUT_BIT_SIZE);
        circuit.output(zero.limbs.clone());

        let dummy_input = create_test_biguint_from_value(ctx.clone(), &params, 0);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &dummy_input, Some(plt_evaluator));

        assert_eq!(eval_result.len(), LIMB_LEN);

        for limb_result in eval_result {
            let limb_coeffs = limb_result.coeffs();
            assert_eq!(*limb_coeffs[0].value(), 0u32.into());
        }
    }

    #[test]
    fn test_biguint_extend_size() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());

        // Extend from 20 bits to 25 bits (5 limbs)
        let extended = big_a.extend_size(25);
        circuit.output(extended.limbs.clone());

        let a_value = create_test_biguint_from_value(ctx.clone(), &params, 12345);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &a_value, Some(plt_evaluator));

        let extended_limb_len = 25 / LIMB_BIT_SIZE + 1; // 6 limbs for 25 bits (with +1 from extend_size)
        assert_eq!(eval_result.len(), extended_limb_len);

        // Check that the original value is preserved in the first limbs
        let mut expected_value = 12345u32;
        let mut expected_limbs = vec![0; extended_limb_len];
        for i in 0..LIMB_LEN {
            if expected_value == 0 {
                break;
            }
            expected_limbs[i] = expected_value % (1u32 << ctx.limb_bit_size);
            expected_value /= 1u32 << ctx.limb_bit_size;
        }
        // Remaining limbs should be 0

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_add_different_limb_sizes() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, LIMB_LEN + 2);

        // Create BigUints with different limb sizes
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b =
            BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..LIMB_LEN + 2].to_vec()); // only 2 limbs
        let result = big_a.add(&big_b, &mut circuit);
        circuit.output(result.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 100);
        let b_limbs = vec![
            DCRTPoly::from_usize_to_constant(&params, 50 % (1u32 << ctx.limb_bit_size) as usize),
            DCRTPoly::from_usize_to_constant(&params, 50 / (1u32 << ctx.limb_bit_size) as usize),
        ];
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b_limbs].concat(),
            Some(plt_evaluator),
        );

        let mut expected_sum = 100u32 + 50u32;
        let mut expected_limbs = vec![0; LIMB_LEN + 1];
        for i in 0..LIMB_LEN + 1 {
            if expected_sum == 0 {
                break;
            }
            expected_limbs[i] = expected_sum % (1u32 << ctx.limb_bit_size);
            expected_sum /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), LIMB_LEN + 1);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_mul_different_limb_sizes() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, LIMB_LEN + 2);

        // Create BigUints with different limb sizes
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b =
            BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..LIMB_LEN + 2].to_vec()); // only 2 limbs
        let result = big_a.mul(&big_b, &mut circuit, None);
        circuit.output(result.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 100);
        let b_limbs = vec![
            DCRTPoly::from_usize_to_constant(&params, 50 % (1u32 << ctx.limb_bit_size) as usize),
            DCRTPoly::from_usize_to_constant(&params, 50 / (1u32 << ctx.limb_bit_size) as usize),
        ];
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b_limbs].concat(),
            Some(plt_evaluator),
        );

        let mut expected_product = 100u32 * 50u32;
        let output_limb_len = (INPUT_BIT_SIZE + 10) / LIMB_BIT_SIZE; // a bit size + b bit size
        let mut expected_limbs = vec![0; output_limb_len];
        for i in 0..output_limb_len {
            if expected_product == 0 {
                break;
            }
            expected_limbs[i] = expected_product % (1u32 << ctx.limb_bit_size);
            expected_product /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), output_limb_len);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_left_shift() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, LIMB_LEN);
        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let shifted = big_a.left_shift(1);
        circuit.output(shifted.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 12345);
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result =
            circuit.eval(&params, &DCRTPoly::const_one(&params), &a, Some(plt_evaluator));

        // Left shift by 1 means removing the first limb
        let mut expected_value = 12345u32;
        let mut expected_limbs = vec![0; LIMB_LEN];
        for i in 0..LIMB_LEN {
            if expected_value == 0 {
                break;
            }
            expected_limbs[i] = expected_value % (1u32 << ctx.limb_bit_size);
            expected_value /= 1u32 << ctx.limb_bit_size;
        }

        // After shifting by 1, we expect to see limbs[1..] from the original
        assert_eq!(eval_result.len(), LIMB_LEN - 1);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i + 1].into());
        }
    }

    #[test]
    fn test_biguint_cmux() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN + 1);

        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b =
            BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..2 * LIMB_LEN].to_vec());
        let selector = inputs[2 * LIMB_LEN];
        let result = big_a.cmux(&big_b, selector, &mut circuit);
        circuit.output(result.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 123);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 456);
        let selector_value = vec![DCRTPoly::from_usize_to_constant(&params, 1)]; // selector = 1, should return 'a'
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b, selector_value].concat(),
            Some(plt_evaluator),
        );

        // With selector = 1, expect to get 'a' (123)
        let mut expected_value = 123u32;
        let mut expected_limbs = vec![0; LIMB_LEN];
        for i in 0..LIMB_LEN {
            if expected_value == 0 {
                break;
            }
            expected_limbs[i] = expected_value % (1u32 << ctx.limb_bit_size);
            expected_value /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), LIMB_LEN);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }

    #[test]
    fn test_biguint_cmux_select_other() {
        let mut circuit = PolyCircuit::<DCRTPoly>::new();
        let (inputs, params, ctx) = create_test_context(&mut circuit, 2 * LIMB_LEN + 1);

        let big_a = BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[0..LIMB_LEN].to_vec());
        let big_b =
            BigUintPoly::<DCRTPoly>::new(ctx.clone(), inputs[LIMB_LEN..2 * LIMB_LEN].to_vec());
        let selector = inputs[2 * LIMB_LEN];
        let result = big_a.cmux(&big_b, selector, &mut circuit);
        circuit.output(result.limbs.clone());

        let a = create_test_biguint_from_value(ctx.clone(), &params, 123);
        let b = create_test_biguint_from_value(ctx.clone(), &params, 456);
        let selector_value = vec![DCRTPoly::from_usize_to_constant(&params, 0)]; // selector = 0, should return 'b'
        let plt_evaluator = PolyPltEvaluator::new();
        let eval_result = circuit.eval(
            &params,
            &DCRTPoly::const_one(&params),
            &[a, b, selector_value].concat(),
            Some(plt_evaluator),
        );

        // With selector = 0, expect to get 'b' (456)
        let mut expected_value = 456u32;
        let mut expected_limbs = vec![0; LIMB_LEN];
        for i in 0..LIMB_LEN {
            if expected_value == 0 {
                break;
            }
            expected_limbs[i] = expected_value % (1u32 << ctx.limb_bit_size);
            expected_value /= 1u32 << ctx.limb_bit_size;
        }

        assert_eq!(eval_result.len(), LIMB_LEN);

        for i in 0..eval_result.len() {
            let coeffs = eval_result[i].coeffs();
            assert_eq!(*coeffs[0].value(), expected_limbs[i].into());
        }
    }
}
