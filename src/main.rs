// helper linear algebra library
extern crate nalgebra;
use nalgebra::{Affine3, Point3, Rotation3, Translation3, Vector2, Vector3};

// add our geo module (exports nalgebra names)
mod geo;
use geo::{RectBounds, Plane};

fn main() {
    // test making a plane
    //let mut plane1 = Plane {
    //    centroid: Point3::origin(),
    //    normal: Vector3::new(0.0, 0.0, 1.0),
    //    bounds: RectBounds{ x: 1.0, y: 2.0 },
    //    global_to_local: Affine3::identity()
    //};

    //let plane2 = Plane {
    //    centroid: Point3::origin(),
    //    normal: Vector3::new(0.0, 0.0, 1.0),
    //    bounds: RectBounds{ x: 1.0, y: 2.0 },
    //    global_to_local: Affine3::identity()
    //};

    //println!("plane centroid starts at: {}", plane1.centroid);

    //// translate plane by -1 and check
    //plane1.translate(Translation3::new(-1.0, 0., 0.));

    //println!("plane centroid ends at: {}", plane1.centroid);

    // make three planes
    let mut pxy = geo::xy_plane( RectBounds::new(1.0, 1.0) );
    let mut pxz = geo::xz_plane( RectBounds::new(1.0, 1.0) );
    let mut pyz = geo::yz_plane( RectBounds::new(1.0, 1.0) );

    //println!("plane1 is {:?}", plane1);
}
