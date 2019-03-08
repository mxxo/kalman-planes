extern crate nalgebra as na;

use na::{Matrix4, Point3, Vector3};

/// The Plane data structure, analagous to PlaneSurface in ACTS
#[derive(Debug)]
struct Plane {
    centroid: Point3<f32>, // global coords
    normal: Vector3<f32>,  // global coords
    bounds: (f32, f32),    // plane coords
    local_to_global: Matrix4<f32> // from the plane to the world
}

fn main() {
    // test making a plane
    let plane1 = Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 0.0, 1.0),
        bounds: (1.0, 2.0),
        local_to_global: Matrix4::identity()
    };

    println!("plane1 is {:?}", plane1);
}
