use num_bigint::BigUint;
use num_traits::{One, Zero};
use rayon::prelude::*;
use std::{env, path::Path};
#[cfg(feature = "disk")]
use tempfile::env::temp_dir;
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
