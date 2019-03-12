extern crate criterion;
extern crate geo;

use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::{Point2, Point3, Rotation3, Translation3, Vector3};

// --

fn bench_stored_inverse(c: &mut Criterion) {
    let pl = dummy_zx_plane();
    c.bench_function("stored inverse", move |b| {
        b.iter(|| pl.get_local_coords(&Point3::new(1002.0, -23.0, 282.)))
    });
}

fn bench_eager_inverse(c: &mut Criterion) {
    let pl = dummy_zx_plane();
    c.bench_function("eager inverse", move |b| {
        b.iter(|| pl.get_local_coords_eager(&Point3::new(1002.0, -23.0, 282.)))
    });
}

// helper for inverse benchmarks
fn dummy_zx_plane() -> geo::Plane {
    geo::zx_plane(geo::RectBounds::new(100.0, 100.0))
        .translate(&Translation3::new(39., -12.0, 2000.))
        .rotate(&Rotation3::from_axis_angle(&Vector3::x_axis(), 0.117117))
        .translate(&Translation3::new(150., 262., -1273.))
}

// --

fn bench_make_xy_plane(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    c.bench_function("make xy plane", move |b| {
        b.iter(|| geo::xy_plane(geo::RectBounds::new(100.0, 100.0)))
    });
}

fn bench_make_zx_plane(c: &mut Criterion) {
    let bounds = geo::RectBounds::new(100.0, 100.0);

    c.bench_function("make zx plane", move |b| {
        b.iter(|| geo::zx_plane(geo::RectBounds::new(100.0, 100.0)))
    });
}

criterion_group!(
    benches,
    bench_stored_inverse,
    bench_eager_inverse,
    bench_make_xy_plane,
    bench_make_zx_plane,
);

criterion_main!(benches);
