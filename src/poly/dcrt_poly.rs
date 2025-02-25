use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPolyImpl},
};
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
    sync::Arc,
};

use super::{params::Params, Polynomial};

use num_traits::Zero;

#[derive(Clone, Debug)]
pub struct DCRTPoly {
    ptr_poly: Arc<UniquePtr<DCRTPolyImpl>>,
}

impl DCRTPoly {
    pub fn new(ptr_poly: UniquePtr<DCRTPolyImpl>) -> Self {
        Self { ptr_poly: ptr_poly.into() }
    }
}

impl Polynomial for DCRTPoly {
    type Error = std::io::Error;

    fn from_const(params: &Params, constant: &u64) -> Result<Self, Self::Error> {
        let res = ffi::DCRTPolyGenFromConst(&params.ptr_params, *constant);
        Ok(DCRTPoly::new(res))
    }
}

// ====== Arithmetic ======

impl Add for DCRTPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let res = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly);
        DCRTPoly::new(res)
    }
}

impl Mul for DCRTPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let res = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly);
        DCRTPoly::new(res)
    }
}

impl Sub for DCRTPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let minus_rhs = rhs.neg();
        self.add(minus_rhs)
    }
}

impl Neg for DCRTPoly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let res = self.ptr_poly.Negate();
        DCRTPoly::new(res)
    }
}

impl PartialEq for DCRTPoly {
    fn eq(&self, _other: &Self) -> bool {
        todo!()
    }
}

impl Eq for DCRTPoly {}

impl Zero for DCRTPoly {
    fn zero() -> Self {
        DCRTPoly::new(UniquePtr::null())
    }

    fn is_zero(&self) -> bool {
        self.ptr_poly.is_null()
    }
}
