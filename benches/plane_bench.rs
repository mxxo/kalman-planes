extern crate criterion;
extern crate geo;

use nalgebra::{Point3, Point2, Translation3, Rotation3, Vector3};
use criterion::{criterion_group, criterion_main, Criterion};

// --

//fn dummy_zx_plane() -> geo::Plane {
//    geo::zx_plane(&geo::RectBounds::new(100.0, 100.0))
//        .translate(&Translation3::new(39., -12.0, 2000))
//        .rotate(&Rotation3::from_axis_angle(
//                &Vector3::x_axis(),
//                0.117117))
//        .translate(&Translation3::new(150, 262, -1273))
//}

fn bench_memoized_inverse(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    let xz_plane = geo::zx_plane(&bounds)
                    .translate(&Translation3::new(39., -12.0, 2000.0))
                    .rotate(&Rotation3::from_axis_angle(
                            &Vector3::x_axis(),
                            0.117117))
                    .translate(&Translation3::new(150., 262., -1273.));

    c.bench_function("memoized inverse", move |b| b.iter(|| xz_plane.get_global_coords(&Point2::new(-100., 200.))));

}

fn bench_eager_inverse(c: &mut Criterion) {}

fn bench_make_xy_plane(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    c.bench_function("make xy plane", move |b| b.iter(|| geo::xy_plane(&bounds)));
}

fn bench_make_zx_plane(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    c.bench_function("make zx plane", move |b| b.iter(|| geo::zx_plane(&bounds)));
}

criterion_group!(
    benches,
    bench_memoized_inverse,
    bench_eager_inverse,
    bench_make_xy_plane,
    bench_make_zx_plane,
);

criterion_main!(benches);
