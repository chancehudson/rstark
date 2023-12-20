use std::rc::Rc;

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use fast_ntt::{
    ntt::{working_modulus, Constants},
    numbers::BigInt as NttBigInt,
    polynomial::Polynomial as NttPolynomial,
};
use itertools::Itertools;
use num_bigint::BigInt;
use rstark::{field::Field, polynomial::Polynomial, BigIntElement, FieldElement};
const DEG: usize = 16;

fn bench_fast_fft_mul<T: FieldElement>(p: &Polynomial<T>, field: &Rc<Field<T>>) {
    Polynomial::fast_mul_fft(p, p, field);
}

fn bench_fft_mul<T: FieldElement>(p: &Polynomial<T>, field: &Rc<Field<T>>) {
    Polynomial::mul_fft(p, p, field);
}

fn bench_mul(x: usize, y: usize, c: &Constants) {
    let ONE = NttBigInt::from(1);
    let a = NttPolynomial::new(vec![0; x].iter().map(|_| ONE).collect_vec());
    let b = NttPolynomial::new(vec![0; y].iter().map(|_| ONE).collect_vec());
    let _ = a.mul(b, c);
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Polynomial Multiplication Benchmarks");

    (2..DEG).for_each(|n| {
        let p = BigIntElement(BigInt::from(3221225473_u32));
        let g = BigIntElement(BigInt::from(5));
        let f = Rc::new(Field::new(p, g));
        let mut poly = Polynomial::new(&f);

        for i in 0..(1 << n) {
            poly.term(&f.random(), i);
        }

        let id = BenchmarkId::new("Old NTT", 1 << n);
        group.bench_function(id, |b| {
            b.iter(|| bench_fft_mul(black_box(&poly), black_box(&f)))
        });

        let id = BenchmarkId::new("New NTT", 1 << n);
        group.bench_function(id, |b| {
            b.iter(|| bench_fast_fft_mul(black_box(&poly), black_box(&f)))
        });

        let id = BenchmarkId::new("Direct NTT", 1 << n);
        let N = NttBigInt::from((2 * n).next_power_of_two());
        let M = N << 1 + 1;
        let c = working_modulus(N, M);
        group.bench_with_input(id, &n, |b, n| {
            b.iter(|| bench_mul(black_box(1 << n), black_box(1 << n), black_box(&c)))
        });
    });
}

criterion_group! {
  name = benches;
  config = Criterion::default().sample_size(10);
  targets = criterion_benchmark
}
criterion_main!(benches);
