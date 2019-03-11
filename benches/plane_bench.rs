extern crate criterion;
use criterion::{Criterion, criterion_group, criterion_main};

// --

fn bench_memoized_inverse(c: &mut Criterion) {
}

fn bench_eager_inverse(c: &mut Criterion) {
}

fn bench_make_xy_plane(c: &mut Criterion) {
    //let xy = geo::Plane::xy_plane(geo::RectBounds::new(2.3, 100.0));
}

fn bench_make_yz_plane(c: &mut Criterion) {

}

criterion_group!(benches, bench_memoized_inverse,
                          bench_eager_inverse,
                          bench_make_xy_plane,
                          bench_make_yz_plane);

criterion_main!(benches);
