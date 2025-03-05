use itertools::Itertools;
use std::ops::{Add, Mul};

use crate::poly::PolyMatrix;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BggPublicKey<M: PolyMatrix> {
    pub matrix: M,
}

impl<M: PolyMatrix> BggPublicKey<M> {
    pub fn new(matrix: M) -> Self {
        Self { matrix }
    }

    pub fn concat_matrix(&self, others: &[Self]) -> M {
        self.matrix.concat_columns(&others.iter().map(|x| x.matrix.clone()).collect_vec()[..])
    }
}

impl<M: PolyMatrix> Add for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        Self { matrix: self.matrix + &other.matrix }
    }
}

impl<M: PolyMatrix> Mul for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggPublicKey<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        let decomposed = other.matrix.decompose();
        let matrix = self.matrix.clone() * decomposed;
        Self { matrix }
    }
}
