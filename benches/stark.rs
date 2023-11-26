#[macro_use]
extern crate criterion;

mod stark_benches {
    use criterion::*;
    use crypto_bigint::U256;
    use num_bigint::BigInt;
    use rstark::{field::Field, BigIntElement, CryptoBigIntElement, FieldElement};

    fn arithmetic(c: &mut Criterion) {
        let mut group: BenchmarkGroup<_> = c.benchmark_group("arithmetic");

        let p = BigIntElement(BigInt::from(1) + BigInt::from(407) * BigInt::from(2).pow(119));
        let g = BigIntElement(BigInt::from(85408008396924667383611388730472331217_u128));

        let f = Field::new(p, g);

        let x = f.bigint(1806576043);
        let y = f.bigint(2045099677);

        let p_2 = CryptoBigIntElement::from_u32(101);
        let g_2 = CryptoBigIntElement::from_u32(0);
        let f_2 = Field::new(p_2, g_2);

        let x_2 = f_2.bigint(1806576043);
        let y_2 = f_2.bigint(2045099677);
        for size in 0..1 {
            group.bench_with_input(
                BenchmarkId::new("add".to_string(), size),
                &size,
                |b, size| {
                    b.iter(|| f.add(&x, &y));
                },
            );

            group.bench_with_input(
                BenchmarkId::new("add with crypto int".to_string(), size),
                &size,
                |b, size| {
                    b.iter(|| f_2.add(&x_2, &y_2));
                },
            );
        }
    }

    criterion_group! {
        name = stark_benches;
        config = Criterion::default();
        targets = arithmetic,
    }
}

criterion_main!(stark_benches::stark_benches,);
