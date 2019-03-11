// helper linear algebra library
extern crate nalgebra;
use nalgebra::{Affine3, Point3, Rotation3, Translation3, Vector2, Vector3};

// add our geo module (exports nalgebra names)
mod geo;
use geo::{RectBounds, Plane};

fn main() {

    //println!("plane centroid starts at: {}", plane1.centroid);

    //// translate plane by -1 and check
    //plane1.translate(Translation3::new(-1.0, 0., 0.));

    //println!("plane centroid ends at: {}", plane1.centroid);

    // make three planes
    let bounds = RectBounds::new(1.0, 1.0);

    let mut pxy = geo::xy_plane( &bounds );
    let mut pxz = geo::xz_plane( &bounds );
    let mut pyz = geo::yz_plane( &bounds );

    let plane_set = vec!(pxy, pxz, pyz);

    //println!("plane1 is {:?}", plane1);
}
