use std::env;
#[cfg(feature = "cpu")]
use std::{thread, time};

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler},
    sampler::{DistType, PolyUniformSampler},
    Poly,
};
use memory_stats::memory_stats;
use num_bigint::BigUint;
use num_traits::{One, Zero};
#[cfg(feature = "cpu")]
use sysinfo::{CpuRefreshKind, RefreshKind, System};
use tracing::{debug, info};

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

pub fn log_mem<T: Into<String>>(tag: T) {
    if let Some(usage) = memory_stats() {
        info!(
            "{} || Current physical/virtual memory usage: {} | {}",
            tag.into(),
            usage.physical_mem,
            usage.virtual_mem,
        );
    } else {
        info!("Couldn't get the current memory usage :(");
    }

    #[cfg(feature = "cpu")]
    {
        let mut sys = System::new_with_specifics(
            RefreshKind::nothing().with_cpu(CpuRefreshKind::everything()),
        );
        sys.refresh_cpu_all();
        thread::sleep(time::Duration::from_millis(200));
        sys.refresh_cpu_all();
        let cpu_usages: Vec<f32> = sys.cpus().iter().map(|cpu| cpu.cpu_usage()).collect();
        info!("CPU usages: {:?} ", cpu_usages);
    }
}

pub fn debug_mem<T: Into<String>>(tag: T) {
    if let Some(usage) = memory_stats() {
        debug!(
            "{} || Current physical/virtual memory usage: {} | {}",
            tag.into(),
            usage.physical_mem,
            usage.virtual_mem
        );
    } else {
        debug!("Couldn't get the current memory usage :(");
    }

    #[cfg(feature = "cpu")]
    {
        let mut sys = System::new_with_specifics(
            RefreshKind::nothing().with_cpu(CpuRefreshKind::everything()),
        );
        sys.refresh_cpu_all();
        thread::sleep(time::Duration::from_millis(200));
        sys.refresh_cpu_all();
        let cpu_usages: Vec<f32> = sys.cpus().iter().map(|cpu| cpu.cpu_usage()).collect();
        debug!("CPU usages: {:?} ", cpu_usages,);
    }
}

pub fn init_tracing() {
    tracing_subscriber::fmt::init();
}

pub fn block_size() -> usize {
    env::var("BLOCK_SIZE").map(|str| str.parse::<usize>().unwrap()).unwrap_or(100)
}

#[macro_export]
macro_rules! parallel_iter {
    ($i: expr) => {{
        rayon::iter::IntoParallelIterator::into_par_iter($i)
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
