use crate::{
    bgg::lut::public_lut::PublicLut,
    poly::{dcrt::DCRTPolyMatrix, Poly, PolyElem, PolyMatrix, PolyParams},
};
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
};

pub trait Evaluable:
    Debug
    + Clone
    + Send
    + Sync
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
{
    type Params: Debug + Clone + Send + Sync;
    type P: Poly;

    fn rotate(self, params: &Self::Params, shift: usize) -> Self;
    fn from_digits(params: &Self::Params, one: &Self, digits: &[u32]) -> Self;
}

impl<P: Poly> Evaluable for P {
    type Params = P::Params;
    type P = P;

    fn rotate(self, params: &Self::Params, shift: usize) -> Self {
        let mut coeffs = self.coeffs();
        coeffs.rotate_right(shift);
        Self::from_coeffs(params, &coeffs)
    }

    fn from_digits(params: &Self::Params, _: &Self, digits: &[u32]) -> Self {
        let coeffs: Vec<P::Elem> = digits
            .par_iter()
            .map(|&digit| <P::Elem as PolyElem>::constant(&params.modulus(), digit as u64))
            .collect();
        Self::from_coeffs(params, &coeffs)
    }
}

pub trait PltEvaluator<E: Evaluable>: Send + Sync {
    fn public_lookup(&self, params: &E::Params, plt: &PublicLut<E::P>, input: E, id: usize) -> E;
}

#[derive(Debug, Clone)]
pub struct PolyPltEvaluator {}
impl<P: Poly> PltEvaluator<P> for PolyPltEvaluator {
    fn public_lookup(&self, _: &P::Params, plt: &PublicLut<P>, input: P, id: usize) -> P {
        // for (_, (x_k, y_k)) in &plt.f {
        //     if *x_k == input {
        //         return y_k.clone();
        //     }
        // }
        // panic!("There is no corresponding input for {id}, {:?}", input.coeffs());
        plt.f[&input].1.clone()
    }
}

impl PolyPltEvaluator {
    pub fn new() -> Self {
        Self {}
    }
}
