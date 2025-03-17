use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler},
    sampler::DistType,
    Poly,
};
use memory_stats::memory_stats;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use tracing::info;

pub fn ceil_log2(q: &BigUint) -> usize {
    assert!(!q.is_zero(), "log2 is undefined for zero");

    let bits = q.bits() as usize;
    if q & (q - BigUint::one()) == BigUint::zero() {
        bits - 1
    } else {
        bits
    }
}

/// Print a ring element
pub fn print_ring_element(label: &str, ring_el: &[u64]) {
    print!("{} [", label);
    for (k, &val) in ring_el.iter().enumerate() {
        if k > 0 {
            print!(", ");
        }
        print!("{}", val);
    }
    println!("]");
}

/// Print a matrix of ring elements
pub fn print_matrix_ring(label: &str, matrix: &[Vec<Vec<u64>>]) {
    println!("\n{}", label,);

    for (i, row) in matrix.iter().enumerate() {
        for (j, col) in row.iter().enumerate() {
            print!("r{}c{}: ", i, j);
            print_ring_element("", col);
        }
    }
}

/// Print a vector of ring elements
pub fn print_vector_ring(label: &str, vec: &[Vec<u64>]) {
    println!("\n{}", label);
    for (i, inner_vec) in vec.iter().enumerate() {
        print!("{}[{}]: ", label, i);
        print_ring_element("", inner_vec);
    }
}

// Helper function to create a random polynomial using UniformSampler
pub fn create_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
    let sampler = DCRTPolyUniformSampler::new();
    sampler.sample_poly(params, &DistType::FinRingDist)
}

pub fn create_bit_random_poly(params: &DCRTPolyParams) -> DCRTPoly {
    let sampler = DCRTPolyUniformSampler::new();
    sampler.sample_poly(params, &DistType::BitDist)
}

// Helper function to create a bit polynomial (0 or 1)
pub fn create_bit_poly(params: &DCRTPolyParams, bit: bool) -> DCRTPoly {
    if bit {
        DCRTPoly::const_one(params)
    } else {
        DCRTPoly::const_zero(params)
    }
}

pub fn log_mem() {
    if let Some(usage) = memory_stats() {
        info!(
            "Current physical/virtural memory usage: {} / {}",
            usage.physical_mem, usage.virtual_mem
        );
    } else {
        info!("Couldn't get the current memory usage :(");
    }
}

pub fn init_tracing() {
    tracing_subscriber::fmt::init();
}

#[macro_export]
macro_rules! parallel_iter {
    ($i: expr) => {{
        #[cfg(not(feature = "parallel"))]
        {
            IntoIterator::into_iter($i)
        }
        #[cfg(feature = "parallel")]
        {
            rayon::iter::IntoParallelIterator::into_par_iter($i)
        }
    }};
}

#[macro_export]
macro_rules! join {
    ($a:expr, $b:expr $(,)?) => {{
        #[cfg(not(feature = "parallel"))]
        {
            ($a(), $b())
        }
        #[cfg(feature = "parallel")]
        {
            rayon::join($a, $b)
        }
    }};
}

/// Implements $tr for all combinations of T and &T by delegating to the &T/&T implementation.
#[macro_export]
macro_rules! impl_binop_with_refs {
    ($T:ty => $tr:ident::$f:ident $($t:tt)*) => {
        impl $tr<$T> for $T {
            type Output = $T;

            #[inline]
            fn $f(self, rhs: $T) -> Self::Output {
                <&$T as $tr<&$T>>::$f(&self, &rhs)
            }
        }

        impl $tr<&$T> for $T {
            type Output = $T;

            #[inline]
            fn $f(self, rhs: &$T) -> Self::Output {
                <&$T as $tr<&$T>>::$f(&self, rhs)
            }
        }

        impl $tr<$T> for &$T {
            type Output = $T;

            #[inline]
            fn $f(self, rhs: $T) -> Self::Output {
                <&$T as $tr<&$T>>::$f(self, &rhs)
            }
        }

        impl $tr<&$T> for &$T {
            type Output = $T;

            #[inline]
            fn $f $($t)*
        }
    };
}
