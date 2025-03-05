use super::BggPublicKey;
use crate::poly::PolyMatrix;
use itertools::Itertools;
use std::ops::{Add, Mul};

#[derive(Debug, Clone)]
pub struct BggEncoding<M: PolyMatrix> {
    pub vector: M,
    pub pubkey: BggPublicKey<M>,
    pub plaintext: Option<<M as PolyMatrix>::P>,
}

impl<M: PolyMatrix> BggEncoding<M> {
    pub fn new(
        vector: M,
        pubkey: BggPublicKey<M>,
        plaintext: Option<<M as PolyMatrix>::P>,
    ) -> Self {
        Self { vector, pubkey, plaintext }
    }

    pub fn concat_vector(&self, others: &[Self]) -> M {
        self.vector.concat_columns(&others.iter().map(|x| x.vector.clone()).collect_vec()[..])
    }
}

impl<M: PolyMatrix> Add for BggEncoding<M> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        self + &other
    }
}

impl<M: PolyMatrix> Add<&Self> for BggEncoding<M> {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        let vector = self.vector + &other.vector;
        let pubkey = self.pubkey + &other.pubkey;
        let plaintext = match (self.plaintext.as_ref(), other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a.clone() + b),
            _ => None,
        };
        Self { vector, pubkey, plaintext }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<M: PolyMatrix> Mul for BggEncoding<M> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<M: PolyMatrix> Mul<&Self> for BggEncoding<M> {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        if self.plaintext.is_none() {
            panic!("Unknown plaintext for the left-hand input of multiplication");
        }
        let decomposed_b = other.pubkey.matrix.decompose();
        let first_term = self.vector.clone() * decomposed_b;
        let second_term = other.vector.clone() * self.plaintext.as_ref().unwrap();
        let new_vector = first_term + second_term;
        let new_plaintext = match (self.plaintext.as_ref(), other.plaintext.as_ref()) {
            (Some(a), Some(b)) => Some(a.clone() * b),
            _ => None,
        };
        let new_pubkey = self.pubkey.clone() * &other.pubkey;
        Self { vector: new_vector, pubkey: new_pubkey, plaintext: new_plaintext }
    }
}
