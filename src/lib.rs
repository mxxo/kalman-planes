/// # A detector geometry module
/// In production, this name would have to be changed to avoid possible
/// name clashes.
///
/// by Max Orok, March 2019

use nalgebra::{Affine3, Point2, Point3, Real, Rotation3, Translation3, Unit, Vector2, Vector3};
use approx::assert_relative_eq;

// EPSILON value for approximate floating point equality
use std::f32;

/// A basic rectangular bounds structure (in local coordinates)
#[derive(Debug, PartialEq)]
pub struct RectBounds{
    pub x: f32,
    pub y: f32,
}

impl RectBounds {
    pub fn new(x: f32, y: f32) -> RectBounds {
        // rough sanity check for zero or negative bounds
        assert!(x > 0.0 && y > 0.0);
        RectBounds { x, y }
    }
}

/// The Plane data structure, somewhat analagous to PlaneSurface in ACTS
#[derive(Debug)]
pub struct Plane<'a> {
    centroid: Point3<f32>, // global coords
    normal: Vector3<f32>,  // global coords
    bounds: &'a RectBounds,    // plane coords (Bounds not owned by Plane)
    local_to_global: Affine3<f32>, // plane to the world
    global_to_local: Option<Affine3<f32>>, // from the world to the plane
                                           // optional because only calculated
                                           // when required by a local_to_global call
}

/// IMPORTANT: only use for debugging!
/// The floating point rough equality checks used here may not be suitable
/// for your application.
///
/// We only compare the geometry and transformation matrix
/// of the two planes i.e. the history of how the plane
/// got to its state doesn't matter

impl<'a> PartialEq for Plane<'a> {
    fn eq(&self, other: &Plane) -> bool {
        self.centroid == other.centroid
        && self.normal.relative_eq(&other.normal, f32::EPSILON, f32::EPSILON)
        && self.bounds == other.bounds // no floating point check since
                                       // bounds can't be scaled in this impl
    }
}

impl<'a> Plane<'a> {

    pub fn get_local_coords(&mut self, global_point: &Point3<f32>) -> Point2<f32> {
        //let local_point_3D = self.global_to_local * global_point;
        // chop off z entry
        //local_point_3D.xy()
        Point2::origin()
    }

    pub fn get_global_coords(&self, local_point : &Point2<f32>) -> Point3<f32> {
        self.local_to_global * Point3::new(local_point.x, local_point.y, 0.0)
    }

    /// expects translation vector in global coords
    /// -> mutates the current local_to_global transform
    pub fn translate(&mut self, translation_vec: &Translation3<f32>) {
        self.centroid = translation_vec * self.centroid;
        self.local_to_global = self.local_to_global * translation_vec ;
    }

    /// expects rotation matrix in global coords
    /// -> mutates the current local_to_global transform
    pub fn rotate(&mut self, rotation_mat: &Rotation3<f32>) {
        self.centroid = rotation_mat * self.centroid;
        self.normal = rotation_mat * self.normal;
        self.local_to_global = rotation_mat * self.local_to_global;
    }
}

/// # Convenience constructors for planes
/// TODO need to think about how the transform matrix should be constructed

/// xy plane has a positive z-direction normal
/// note: local to global transform is an identity matrix (l_x = g_x; l_y = g_y)
pub fn xy_plane(bounds: &RectBounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 0.0, 1.0),
        bounds,
        local_to_global: Affine3::identity(),
        global_to_local: None,
    }
}

/// zx plane has a positive y-direction normal
/// note: reuse xy cstor and make a zx plane by
/// rotating the xy plane about the x-axis by -pi / 2 radians
pub fn zx_plane(bounds: &RectBounds) -> Plane {

    let mut pl = xy_plane(&bounds);
    pl.rotate(
        &Rotation3::from_axis_angle(&Vector3::x_axis(),
        -1.0 * f32::consts::FRAC_PI_2)
    );

    pl
}

/// yz plane has a positive x-direction normal
/// note: reuse xy cstor and make a yz plane by
/// rotating the xy plane about the y-axis
pub fn yz_plane(bounds: &RectBounds) -> Plane {

    let mut pl = xy_plane(&bounds);
    pl.rotate(
        &Rotation3::from_axis_angle(&Vector3::y_axis(),
        Real::frac_pi_2())
    );

    pl
}



/// Unit tests for planes

#[cfg(test)]
mod tests {
    // import all names from the outer scope
    use super::*;

    /// RectBounds tests
    #[test]
    #[should_panic]
    fn positive_bounds() {
       RectBounds::new(-1.0, 0.0);
    }

    #[test]
    #[should_panic]
    fn zero_bounds() {
        RectBounds::new(0.0, 1.0);
    }

    /// Local and global coordinate systems

    #[test]
    fn local_to_global() {

        let bounds = RectBounds::new(1.0, 1.0);
        let mut pxy = xy_plane( &bounds );

        let translation_vec = Vector3::new(1.0, -1.0, 3.0);
        let global_point = Point3::from(translation_vec);

        pxy.translate(&Translation3::from(translation_vec));

        // point at the origin of the plane should be at the
        // translation point
        let converted_global = pxy.get_global_coords(&Point2::origin());
        assert_relative_eq!(global_point, converted_global);
    }

    #[test]
    /// tests coordinate conversion after a translation and rotation
    fn global_to_local() {
        //let bounds = RectBounds::new(1.0, 1.0);

        //let mut pxy = xy_plane( &bounds );
        //let translation_vec = Vector3::new(1.0, -1.0, 3.0);
        //let global_point = Point3::from(translation_vec);

        //pxy.translate(&Translation3::from(translation_vec));
        //pxy.rotate(&Rotation3::from_axis_angle(&Vector3::y_axis(), Real::frac_pi_2()));
        //// point at the origin of the plane should be at the
        //// translation point
        //let local_point = pxy.get_local_coords(&global_point);
        //assert_relative_eq!(local_point, Point2::origin());
        assert!(false);
    }

    #[test]
    fn xy_plane_coords() {

        let bounds = RectBounds::new(1.0, 1.0);
        let pl_xy = xy_plane(&bounds);

        // test x-axis points
        let globl_x = Point3::new(1.0, 0.0, 0.0);
        let local_x = Point2::new(1.0, 0.0);

        assert_relative_eq!(globl_x,
                            pl_xy.get_global_coords(&local_x));

        // test y-axis points
        let globl_y = Point3::new(0.0, 1.0, 0.0);
        let local_y = Point2::new(0.0, 1.0);

        assert_relative_eq!(globl_y,
                            pl_xy.get_global_coords(&local_y));

        // test x and y
        let globl_xy = Point3::new(-3.0, -3.0, 0.0);
        let local_xy = Point2::new(-3.0, -3.0);

        assert_relative_eq!(globl_xy,
                            pl_xy.get_global_coords(&local_xy));
    }

    #[test]
    fn zx_plane_coords() {
        let bounds = RectBounds::new(1.0, 1.0);
        let pl_zx = zx_plane(&bounds);

        // y-points should be ignored
//        let y_point = Point3::new(0.0, 1.0, 0.0);
//        assert_relative_eq!(y_point,
//                            pl_zx.get_global_coords(&Point2::new(0.0, 0.0)));
//

        // z-axis becomes local x-axis -> z point should equal local x point
        let globl_z = Point3::new(1.0, 0.0, 1.0);
        let local_z = Point2::new(1.0, 1.0);

        assert_relative_eq!(globl_z,
                            pl_zx.get_global_coords(&local_z));
        // on the xz plane, local x is global x, and local y is global z
        // therefore, local: (2.0, -1.0) => global: (2.0, 0.0, -1.0)

        //let global_point = Point3::new(2.0, 0.0, -1.0);

        //assert_relative_eq!(global_point,
        //                    pl_zx.get_global_coords(&Point2::new(2.0, -1.0)));
    }

    #[test]
    fn yz_plane_coords() {
        assert!(false);
        let bounds = RectBounds::new(1.0, 1.0);
    }

    /// Transformation tests
    #[test]
    fn reversible_translation() {

        let bounds = RectBounds::new(1.0, 2.0);

        let mut p1 = xy_plane( &bounds );
        let p2 = xy_plane( &bounds );

        // translate p1 there and back and test equality
        assert_eq!(p1, p2);

        p1.translate(&Translation3::new(1.0, 2.0, 3.0));
        p1.translate(&Translation3::new(-1.0, -2.0, -3.0));

        assert_eq!(p1, p2);
    }

    #[test]
    fn reversible_rotation() {

        let bounds = RectBounds::new(3.0, 3.0);

        // two planes centred at the origin
        let mut p1 = xy_plane ( &bounds );
        let p2 = yz_plane ( &bounds );

        // rotate the xy_plane pi/2 radians around y
        // to become identical to a xz_plane
        let y_axis: Unit<Vector3<f32>> = Vector3::y_axis();
        let y_rot = Rotation3::from_axis_angle(&y_axis, Real::frac_pi_2());

        p1.rotate(&y_rot);

        // note: approx equal for floating point numbers!
        assert_eq!(p1, p2);
    }

    // test that global from local and back again gives the same vector
    #[test]
    fn reversible_coord_transform() {
        assert_eq!(1, 2);
    }
}
