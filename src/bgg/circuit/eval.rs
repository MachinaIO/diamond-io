use crate::{
    bgg::lut::public_lut::PublicLut,
    poly::{dcrt::DCRTPolyMatrix, Poly, PolyElem, PolyMatrix, PolyParams},
};
use rayon::prelude::*;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Sub},
    path::PathBuf,
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
    type Matrix: PolyMatrix;

    fn rotate(self, params: &Self::Params, shift: usize) -> Self;
    fn from_digits(params: &Self::Params, one: &Self, digits: &[u32]) -> Self;
    fn public_lookup(
        self,
        params: &Self::Params,
        plt: &mut PublicLut<Self::Matrix>,
        helper: Option<(Self::Matrix, PathBuf, usize, usize)>,
    ) -> Self;
}

impl<P: Poly> Evaluable for P {
    type Params = P::Params;
    type Matrix = DCRTPolyMatrix;

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

    fn public_lookup(
        self,
        _: &Self::Params,
        _: &mut PublicLut<Self::Matrix>,
        _: Option<(Self::Matrix, PathBuf, usize, usize)>,
    ) -> Self {
        self
    }
}
