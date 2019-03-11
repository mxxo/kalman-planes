extern crate nalgebra;  // helper linear algebra library
extern crate approx;    // relative_eq! macros for floating point comparisons

use nalgebra::{Affine3, Point2, Point3, Rotation3, Translation3, Vector2, Vector3};

// add our geo module
extern crate geo;
use geo::{RectBounds, Plane};

fn main() {

    //println!("plane centroid starts at: {}", plane1.centroid);

    //// translate plane by -1 and check
    //plane1.translate(Translation3::new(-1.0, 0., 0.));

    //println!("plane centroid ends at: {}", plane1.centroid);

    // make three planes
    let bounds = RectBounds::new(1.0, 1.0);

    let mut pxy = geo::xy_plane( &bounds );

    pxy.translate(&Translation3::new(1.0, -1.0, 3.0));

    //let glob = Point3::new(1.0, 2.0, 3.0);
    //pxy.get_local_coords(&glob);

    //println!("local coords are {:?}", pxy.get_local_coords(&glob));
    //println!("global coords are {:?}", glob);

    let loc = Point2::origin();
    println!("local coords are {:?}", &loc);
    println!("global coords are {:?}", pxy.get_global_coords(&loc));

    println!("plane is {:?}", pxy);

    let mut pxz = geo::xz_plane( &bounds );
    let mut pyz = geo::yz_plane( &bounds );



    // let plane_set = vec!(pxy, pxz, pyz);

    //println!("plane1 is {:?}", plane1);
}
