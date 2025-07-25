use crate::{
    bgg::{
        circuit::{GateId, PolyCircuit},
        lut::public_lut::PublicLut,
    },
    poly::Poly,
};
use std::{collections::HashMap, marker::PhantomData, sync::Arc};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BigUintContext<P: Poly> {
    pub limb_bit_size: usize,
    pub const_zero: GateId,
    pub const_base: GateId,
    pub add_lut_ids: (usize, usize),
    pub mul_lut_ids: (usize, usize),
    _p: PhantomData<P>,
}

impl<P: Poly> BigUintContext<P> {
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
            let input = P::const_int(params, k);
            let output_f = P::const_int(params, k % base);
            let output_g = P::const_int(params, k / base);
            f.insert(input.clone(), (k, output_f));
            g.insert(input, (k, output_g));
        }
        (PublicLut::new(f), PublicLut::new(g))
    }
}

#[derive(Debug, Clone)]
pub struct BigUint<P: Poly> {
    pub ctx: Arc<BigUintContext<P>>,
    pub limbs: Vec<GateId>,
    _p: PhantomData<P>,
}

impl<P: Poly> BigUint<P> {
    pub fn new(ctx: Arc<BigUintContext<P>>, limbs: Vec<GateId>) -> Self {
        Self { ctx, limbs, _p: PhantomData }
    }

    pub fn bit_size(&self) -> usize {
        self.limbs.len() * self.ctx.limb_bit_size
    }

    pub fn zero(ctx: Arc<BigUintContext<P>>, bit_size: usize) -> Self {
        debug_assert_eq!(bit_size % ctx.limb_bit_size, 0);
        let limb_len = bit_size / ctx.limb_bit_size;
        let limbs = vec![ctx.const_zero; limb_len];
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

    pub fn add(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> Self {
        debug_assert_eq!(self.limbs.len(), other.limbs.len());
        debug_assert_eq!(self.ctx, other.ctx);
        let mut limbs = Vec::with_capacity(self.limbs.len() + 1);
        let mut carry = circuit.const_zero_gate();
        for i in 0..self.limbs.len() {
            let tmp = circuit.add_gate(self.limbs[i], other.limbs[i]);
            let sum = circuit.add_gate(tmp, carry);
            let sum_l = circuit.public_lookup_gate(sum, self.ctx.add_lut_ids.0);
            let sum_h = circuit.public_lookup_gate(sum, self.ctx.add_lut_ids.1);
            limbs.push(sum_l);
            carry = sum_h;
        }
        limbs.push(carry);
        Self { ctx: self.ctx.clone(), limbs, _p: PhantomData }
    }

    pub fn geq(&self, other: &Self, circuit: &mut PolyCircuit<P>) -> (GateId, Self) {
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
        debug_assert_eq!(self.limbs.len(), other.limbs.len());
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
}
