use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use exact_number::{based_expr, rat::{Integer, Rat}, BasedExpr};

pub fn bench_add_small_integers(c: &mut Criterion) {
    let data = [based_expr!(12345678), based_expr!(87654321)];
    let baseless_data = [BasedExpr::Baseless(12345678.into()), BasedExpr::Baseless(87654321.into())];
    let rat_data = [Rat::from(12345678), Rat::from(87654321)];
    let int_data = [Integer::from(12345678), Integer::from(87654321)];
    c.bench_function("bench_add_small_integers_int", move |b| {
        b.iter_batched(|| int_data.clone(), |[a, b]| a + b, BatchSize::SmallInput);
    });
    c.bench_function("bench_add_small_integers_rat", move |b| {
        b.iter_batched(|| rat_data.clone(), |[a, b]| a + b, BatchSize::SmallInput);
    });
    c.bench_function("bench_add_small_baseless", move |b| {
        b.iter_batched(|| baseless_data.clone(), |[a, b]| a + b, BatchSize::SmallInput);
    });
    c.bench_function("bench_add_small_integers", move |b| {
        b.iter_batched(|| data.clone(), |[a, b]| a + b, BatchSize::SmallInput);
    });
}

criterion_group!(add, bench_add_small_integers);
criterion_main!(add);