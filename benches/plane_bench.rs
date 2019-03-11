extern crate criterion;
use criterion::{Criterion, criterion_group, criterion_main};

extern crate geo;

// --

fn bench_memoized_inverse(c: &mut Criterion) {
    //c.bench_function("inverse global_to_local ", |b| b.iter(||
}

fn bench_eager_inverse(c: &mut Criterion) {

}

fn bench_make_xy_plane(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    c.bench_function("make xy plane",
        move |b| b.iter(|| geo::xy_plane(&bounds))
    );
}

fn bench_make_zx_plane(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    c.bench_function("make zx plane",
        move |b| b.iter(|| geo::zx_plane(&bounds))
    );
}

criterion_group!(benches, bench_memoized_inverse,
                          bench_eager_inverse,
                          bench_make_xy_plane,
                          bench_make_zx_plane,
                );

criterion_main!(benches);
