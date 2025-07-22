use std::{env, path::Path};
#[cfg(feature = "cpu")]
use std::{thread, time};

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler},
    sampler::{DistType, PolyUniformSampler},
    Poly, PolyMatrix,
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
use tokio::runtime::Handle;
use tracing::{debug, info};
use walkdir::WalkDir;

pub fn store_and_drop_matrix<M>(matrix: M, dir: &Path, id: &str) -> tokio::task::JoinHandle<()>
where
    M: PolyMatrix + Send + 'static,
{
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    tokio::task::spawn_blocking(move || {
        log_mem(format!("Storing {id}"));
        Handle::current().block_on(async {
            matrix.write_to_files(&dir, &id).await;
        });
        drop(matrix);
        log_mem(format!("Stored {id}"));
    })
}

#[cfg(feature = "debug")]
pub fn store_and_drop_poly<P: Poly + 'static>(
    poly: P,
    dir: &Path,
    id: &str,
) -> tokio::task::JoinHandle<()> {
    let dir = dir.to_path_buf();
    let id = id.to_owned();

    tokio::spawn(async move {
        log_mem(format!("Storing {id}"));
        poly.write_to_file(&dir, &id).await;
        drop(poly);
        log_mem(format!("Stored {id}"));
    })
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

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::poly::dcrt::DCRTPolyMatrix;
//     use futures::future::join_all;
//     use std::{fs, path::Path, time::Instant};
//     use tokio::{sync::mpsc, task::JoinSet};
//     use tracing::info;

//     /*
//          ┌──────────────── comparison ───────────────┐
//     2025-07-22T06:44:05.495608Z  INFO diamond_io::utils::tests: │ vec + join_all  : 114426 ms │
//     2025-07-22T06:44:05.495630Z  INFO diamond_io::utils::tests: │ pipeline mpsc   : 109998 ms │
// (faster 1.0×)     2025-07-22T06:44:05.495635Z  INFO diamond_io::utils::tests:
// └───────────────────────────────────────────── */     /// Adjust to taste
//     const N_MATS: usize = 20; // how many matrices we write
//     const CHANNEL_CAP: usize = 5; // queue depth for the pipeline

//     #[tokio::test(flavor = "multi_thread", worker_threads = 4)]
//     async fn compare_vec_vs_pipeline() {
//         init_tracing();
//         std::env::set_var("BLOCK_SIZE", "1000");

//         // ─── crypto params & sampler ───────────────────────────────────────────────
//         let params = DCRTPolyParams::new(8192, 7, 51, 17);
//         let sampler = DCRTPolyUniformSampler::new();
//         let dist = DistType::BitDist;

//         // ─── generate the matrices *once* so both variants write the same data ——──
//         let mats: Vec<_> = (0..N_MATS)
//             .map(|i| (format!("mat_{i}"), sampler.sample_uniform(&params, 32, 32, dist)))
//             .collect();

//         // fresh output dirs
//         let dir_vec = Path::new("out_vec");
//         let dir_pipe = Path::new("out_pipe");
//         for d in [dir_vec, dir_pipe] {
//             if d.exists() {
//                 fs::remove_dir_all(d).unwrap();
//             }
//             fs::create_dir(d).unwrap();
//         }

//         // ╔════════════════════════════════════════════╗
//         // ║ 1. current vec + join_all                  ║
//         // ╚════════════════════════════════════════════╝
//         let t0_vec = Instant::now();
//         let mut handles = Vec::new();

//         for (id, mat) in mats.iter().cloned() {
//             handles.push(store_and_drop_matrix(mat, dir_vec, &id));
//         }

//         // wait for them all
//         join_all(handles).await;
//         let dur_vec = t0_vec.elapsed();

//         // ╔════════════════════════════════════════════╗
//         // ║ 2. pipeline with bounded channel           ║
//         // ╚════════════════════════════════════════════╝
//         let t0_pipe = Instant::now();

//         let (tx, mut rx) = mpsc::channel::<(String, DCRTPolyMatrix)>(CHANNEL_CAP);

//         // single writer task; uses JoinSet so matrix writes overlap
//         let dir_pipe_clone = dir_pipe.to_path_buf();
//         let writer = tokio::spawn(async move {
//             let mut js = JoinSet::new();
//             while let Some((id, mat)) = rx.recv().await {
//                 js.spawn(store_and_drop_matrix(mat, &dir_pipe_clone, &id));
//             }
//             while let Some(r) = js.join_next().await {
//                 r.expect("matrix‑store failed").unwrap();
//             }
//         });

//         // producer side
//         for (id, mat) in mats.iter().cloned() {
//             tx.send((id, mat)).await.unwrap();
//         }
//         drop(tx); // close channel
//         writer.await.unwrap();

//         let dur_pipe = t0_pipe.elapsed();

//         // ─── perf summary ──────────────────────────────────────────────────────────
//         info!("┌──────────────── comparison ───────────────┐");
//         info!("│ vec + join_all  : {:>6} ms │", dur_vec.as_millis());
//         info!(
//             "│ pipeline mpsc   : {:>6} ms │ ({}{:.1}×)",
//             dur_pipe.as_millis(),
//             if dur_pipe < dur_vec { "faster " } else { "slower " },
//             dur_vec.as_secs_f64() / dur_pipe.as_secs_f64()
//         );
//         info!("└─────────────────────────────────────────────┘");

//         // sanity: pipeline should not be catastrophically slower
//         assert!(dur_pipe < dur_vec * 2, "pipeline unexpectedly slower than vec/join_all");
//     }
// }
