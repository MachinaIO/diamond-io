use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use diamond_io::poly::dcrt::{DCRTPolyMatrix, DCRTPolyParams};

fn bench_matrix_operation(c: &mut Criterion) {
    let params = DCRTPolyParams::new(4, 2, 17, 1);
    let size = 3;
    let matrix_a = DCRTPolyMatrix::identity(&params, size, None);
    let matrix_b = DCRTPolyMatrix::identity(&params, size, None);

    c.bench_with_input(BenchmarkId::new("Matrix Addition", size), &size, |b, _| {
        b.iter(|| {
            let _ = matrix_a.clone() + &matrix_b;
        })
    });

    c.bench_with_input(BenchmarkId::new("Matrix Multiplication", size), &size, |b, _| {
        b.iter(|| {
            let _ = matrix_a.clone() * &matrix_b;
        })
    });
}

fn bench_matrix_tensor(c: &mut Criterion) {
    let params = DCRTPolyParams::new(4, 2, 17, 1);
    let shape_a = (2, 2);
    let shape_b = (3, 3);

    let matrix_a = DCRTPolyMatrix::identity(&params, 2, None);
    let matrix_b = DCRTPolyMatrix::identity(&params, 3, None);

    c.bench_with_input(
        BenchmarkId::new("Matrix Tensor", format!("{shape_a:?} ⊗ {shape_b:?}")),
        &(matrix_a, matrix_b),
        |b, (m_a, m_b)| {
            b.iter(|| {
                let _ = m_a.tensor(m_b);
            })
        },
    );
}

criterion_group!(benches, bench_matrix_operation, bench_matrix_tensor);
criterion_main!(benches);
