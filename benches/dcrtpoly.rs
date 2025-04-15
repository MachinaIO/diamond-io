use criterion::{criterion_group, criterion_main, Criterion};
use diamond_io::poly::{
    dcrt::{DCRTPoly, DCRTPolyParams, FinRingElem},
    Poly, PolyParams,
};

pub fn dcrtpoly_operations_bench(c: &mut Criterion) {
    let params = DCRTPolyParams::new(4, 2, 17, 1);
    let q = params.modulus();

    let coeffs1: Vec<FinRingElem> = (0..params.ring_dimension())
        .map(|_| FinRingElem::new(rand::random::<u32>(), q.clone()))
        .collect();
    let coeffs2: Vec<FinRingElem> = (0..params.ring_dimension())
        .map(|_| FinRingElem::new(rand::random::<u32>(), q.clone()))
        .collect();

    let poly1 = DCRTPoly::from_coeffs(&params, &coeffs1);
    let poly2 = DCRTPoly::from_coeffs(&params, &coeffs2);

    c.bench_function("DCRTPoly add", |b| {
        b.iter(|| {
            let _result = &poly1 + &poly2;
        })
    });

    c.bench_function("DCRTPoly mul", |b| {
        b.iter(|| {
            let _result = &poly1 * &poly2;
        })
    });

    c.bench_function("DCRTPoly sub", |b| {
        b.iter(|| {
            let _result = &poly1 - &poly2;
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = dcrtpoly_operations_bench
);
criterion_main!(benches);
