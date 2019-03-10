/// # A detector geometry module
/// In production, this name would have to be changed to avoid possible
/// name clashes.
/// by Max Orok, March 2019

use nalgebra::{Affine3, Point3, Real, Rotation3, Translation3, Unit, Vector2, Vector3};

// EPSILON values
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
pub struct Plane {
    pub centroid: Point3<f32>, // global coords
    pub normal: Vector3<f32>,  // global coords
    pub bounds: RectBounds,    // plane coords (Bounds object is owned by Plane)
    pub global_to_local: Affine3<f32> // from the world to the plane
}


/// IMPORTANT: only use for debugging!
/// The floating point rough equality checks used here may not be suitable
/// for your application.
///
/// We only compare the geometry and transformation matrix
/// of the two planes i.e. the history of how the plane
/// got to its state doesn't matter

impl PartialEq for Plane {
    fn eq(&self, other: &Plane) -> bool {
        self.centroid == other.centroid
        && self.normal.relative_eq(&other.normal, f32::EPSILON, f32::EPSILON)
        && self.bounds == other.bounds // no floating point check since
                                       // bounds can't be scaled in this impl

        // comment out for now since it might not be use and
        // floating point errors can mess this up
        // && self.global_to_local == other.global_to_local
    }
}

impl Plane {

    pub fn get_local_coords(&self, global_vec: &Vector3<f32>) -> Vector2<f32> {
        // broken impl for now
        Vector2::new(0.0, 0.0)
    }

    /// expects translation vector in global coords
    /// we choose to consume the translation here
    /// -> mutates the current global_to_local transform
    pub fn translate(&mut self, translation_vec: &Translation3<f32>) {
        self.centroid = translation_vec * self.centroid;
        self.global_to_local = translation_vec * self.global_to_local;
    }

    /// expects rotation matrix in global coords
    /// -> mutates the current global_to_local transform
    pub fn rotate(&mut self, rotation_mat : &Rotation3<f32>) {
        self.centroid = rotation_mat * self.centroid;
        self.normal = rotation_mat * self.normal;
        self.global_to_local = rotation_mat * self.global_to_local;
    }
}

/// # Convenience constructors for planes
/// TODO need to think about how the transform matrix should be constructed

/// xy plane has a positive z-direction normal
pub fn xy_plane(bounds: RectBounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 0.0, 1.0),
        bounds,
        global_to_local: Affine3::identity()
    }
}

/// xz plane has a positive y-direction normal
pub fn xz_plane(bounds: RectBounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 1.0, 0.0),
        bounds,
        global_to_local: Affine3::identity()
    }
}

/// yz plane has a positive x-direction normal
pub fn yz_plane(bounds: RectBounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(1.0, 0.0, 0.0),
        bounds,
        global_to_local: Affine3::identity()
    }
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

    /// Transformation tests
    #[test]
    fn reversible_translation() {

        let mut p1 = xy_plane( RectBounds::new(1.0, 2.0) );
        let p2 = xy_plane( RectBounds::new(1.0, 2.0) );

        // translate p1 there and back and test equality
        assert_eq!(p1, p2);

        p1.translate(&Translation3::new(1.0, 2.0, 3.0));
        p1.translate(&Translation3::new(-1.0, -2.0, -3.0));

        assert_eq!(p1, p2);
    }

    #[test]
    fn reversible_rotation() {

        // two planes centred at the origin
        let mut p1 = xy_plane ( RectBounds::new(3.0, 3.0) );
        let p2 = yz_plane ( RectBounds::new(3.0, 3.0) );

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
