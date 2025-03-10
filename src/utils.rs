use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, DCRTPolyUniformSampler},
    sampler::DistType,
};

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
