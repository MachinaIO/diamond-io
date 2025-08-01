use std::{env, path::Path};
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
use rayon::prelude::*;
use std::time::{Duration, Instant};
#[cfg(feature = "cpu")]
use sysinfo::{CpuRefreshKind, RefreshKind, System};
#[cfg(feature = "disk")]
use tempfile::env::temp_dir;
use tracing::{debug, info};
use walkdir::WalkDir;

/// ideal thread chunk size for parallel
pub fn chunk_size_for(original: usize) -> usize {
    match original {
        0..=2048 => 128,
        2049..=4096 => 256,
        4097..=8192 => 512,
        _ => 1024,
    }
}

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
    print!("{label} [");
    for (k, &val) in ring_el.iter().enumerate() {
        if k > 0 {
            print!(", ");
        }
        print!("{val}");
    }
    println!("]");
}

/// Print a matrix of ring elements
pub fn print_matrix_ring(label: &str, matrix: &[Vec<Vec<u64>>]) {
    println!("\n{label}",);

    for (i, row) in matrix.iter().enumerate() {
        for (j, col) in row.iter().enumerate() {
            print!("r{i}c{j}: ");
            print_ring_element("", col);
        }
    }
}

/// Print a vector of ring elements
pub fn print_vector_ring(label: &str, vec: &[Vec<u64>]) {
    println!("\n{label}");
    for (i, inner_vec) in vec.iter().enumerate() {
        print!("{label}[{i}]: ");
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

pub fn timed_read<T, F: FnOnce() -> T>(label: &str, f: F, total: &mut Duration) -> T {
    let start = Instant::now();
    let res = f();
    let elapsed = start.elapsed();
    *total += elapsed;
    crate::utils::log_mem(format!("{label} loaded in {elapsed:?}"));
    res
}

/// Execute a function and log the elapsed time.
pub fn timed<T, F: FnOnce() -> T>(label: &str, f: F) -> T {
    let start = Instant::now();
    let res = f();
    let elapsed = start.elapsed();
    info!("{} completed in {:?}", label, elapsed);
    res
}

pub fn init_tracing() {
    tracing_subscriber::fmt::init();
}

pub fn block_size() -> usize {
    env::var("BLOCK_SIZE").map(|str| str.parse::<usize>().unwrap()).unwrap_or(100)
}

/// Calculate the total size of a directory in bytes
pub fn calculate_directory_size<P: AsRef<Path>>(path: P) -> u64 {
    WalkDir::new(path)
        .follow_links(false)
        .into_iter()
        .par_bridge()
        .filter_map(Result::ok)
        .filter_map(|e| e.metadata().ok())
        .filter(|m| m.is_file())
        .map(|m| m.len())
        .sum()
}

#[cfg(feature = "disk")]
pub fn calculate_tmp_size() -> u64 {
    calculate_directory_size(temp_dir())
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
