use criterion::{black_box, criterion_group, criterion_main, Criterion};
use diamond_io::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, FinRingElem},
    Poly, PolyParams,
};
use num_bigint::BigUint;

fn bench_decompose_dim_32(c: &mut Criterion) {
    let params = DCRTPolyParams::new(32, 16, 60);
    let q = params.modulus();

    let poly1 = DCRTPoly::from_coeffs(
        &params,
        &vec![
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(400u32), q.clone()),
        ],
    );
    let poly2 = poly1.decompose(&params);

    c.bench_function("DCRTPoly decompose dim 32", |b| b.iter(|| black_box(poly2.clone())));
}

fn bench_addition_dim_4(c: &mut Criterion) {
    let params = DCRTPolyParams::new(4, 2, 60);
    let q = params.modulus();

    let poly1 = DCRTPoly::from_coeffs(
        &params,
        &[
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
            FinRingElem::new(BigUint::from(300u32), q.clone()),
        ],
    );

    c.bench_function("DCRTPoly decompose dim 4", |b| b.iter(|| black_box(poly1.clone())));
}

fn bench_multiplication(c: &mut Criterion) {
    let params = DCRTPolyParams::new(4, 2, 17);
    let q = params.modulus();

    let poly1 = DCRTPoly::from_coeffs(
        &params,
        &[
            FinRingElem::new(BigUint::from(100u32), q.clone()),
            FinRingElem::new(BigUint::from(200u32), q.clone()),
        ],
    );
    poly1.decompose(&params);

    c.bench_function("DCRTPoly multiplication", |b| b.iter(|| black_box(poly1.clone())));
}

fn bench_decomposition(c: &mut Criterion) {
    let params = DCRTPolyParams::new(4, 2, 17);
    let poly = DCRTPoly::const_power_of_two(&params, 5);

    c.bench_function("DCRTPoly decomposition", |b| b.iter(|| black_box(poly.decompose(&params))));
}

criterion_group!(
    benches,
    bench_addition_dim_4,
    bench_decompose_dim_32,
    bench_multiplication,
    bench_decomposition
);
criterion_main!(benches);
