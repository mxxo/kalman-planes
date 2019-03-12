/// # A detector geometry module
/// In production, this name would have to be changed to avoid possible
/// name clashes.
///
/// by Max Orok, March 2019
use approx::assert_relative_eq;
use nalgebra::{Affine3, Point2, Point3, Real, Rotation3, Translation3, Unit, Vector2, Vector3};

// EPSILON value for approximate floating point equality
use std::f32;

/// A basic rectangular bounds structure (in local coordinates)
#[derive(Debug, PartialEq)]
pub struct RectBounds {
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
    centroid: Point3<f32>,         // global coords
    normal: Vector3<f32>,          // global coords
    bounds: RectBounds,            // plane coords (owned by Plane)
    local_to_global: Affine3<f32>, // plane to the world
    global_to_local: Affine3<f32>, // from the world to the plane
                                   // more expensive than local_to_global
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
            && self
                .normal
                .relative_eq(&other.normal, f32::EPSILON, f32::EPSILON)
            && self.bounds == other.bounds // no floating point check since
                                           // bounds can't be scaled in this impl
    }
}

impl Plane {
    // takes mut ref to self because
    pub fn get_local_coords(&self, global_point: &Point3<f32>) -> Point2<f32> {
        //let local_point_3D = match self.global_to_local {
        //    Some(inverse_matrix) => inverse_matrix * global_point,
        //    None => {
        //        let inv = self.local_to_global.inverse();
        //        self.global_to_local = Some(inv);
        //        inv * global_point
        //    }
        //};
        let local_point3D = self.global_to_local * global_point;

        // chop off z entry
        local_point3D.xy()
    }

    pub fn get_local_coords_eager(&self, global_point: &Point3<f32>) -> Point2<f32> {
        let local_point_3D = self.local_to_global.inverse() * global_point;
        local_point_3D.xy()
    }

    pub fn get_global_coords(&self, local_point: &Point2<f32>) -> Point3<f32> {
        self.local_to_global * Point3::new(local_point.x, local_point.y, 0.0)
    }

    /// expects translation vector in global coords
    /// -> mutates the current local_to_global transform
    pub fn translate(mut self, translation_vec: &Translation3<f32>) -> Self {
        self.centroid = translation_vec * self.centroid;
        self.local_to_global = self.local_to_global * translation_vec;
        // update current inverse matrix
        self.global_to_local = self.local_to_global.inverse();

        self
    }

    /// expects rotation matrix in global coords
    /// -> mutates the current local_to_global transform
    pub fn rotate(mut self, rotation_mat: &Rotation3<f32>) -> Self {
        self.centroid = rotation_mat * self.centroid;
        self.normal = rotation_mat * self.normal;
        self.local_to_global = rotation_mat * self.local_to_global;
        // nullify current inverse matrix
        self.global_to_local = self.local_to_global.inverse();

        self
    }
}

/// # Convenience constructors for planes
/// TODO need to think about how the transform matrix should be constructed

/// xy plane has a positive z-direction normal
/// note: local to global transform is an identity matrix (l_x = g_x; l_y = g_y)
pub fn xy_plane(bounds: RectBounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 0.0, 1.0),
        bounds,
        local_to_global: Affine3::identity(),
        global_to_local: Affine3::identity(),
    }
}

/// zx plane has a positive y-direction normal
/// note: reuse xy cstor and make a zx plane by
/// rotating the xy plane about the x-axis by -pi / 2 radians
pub fn zx_plane(bounds: RectBounds) -> Plane {
    xy_plane(bounds)
        .rotate(&Rotation3::from_axis_angle(
            &Vector3::x_axis(),
            -1.0 * f32::consts::FRAC_PI_2,
        ))
        .rotate(&Rotation3::from_axis_angle(
            &Vector3::y_axis(),
            -1.0 * f32::consts::FRAC_PI_2,
        ))
}

/// yz plane has a positive x-direction normal
/// note: reuse xy cstor and make a yz plane by
/// rotating the xy plane about the y-axis
pub fn yz_plane(bounds: RectBounds) -> Plane {
    xy_plane(bounds)
        .rotate(&Rotation3::from_axis_angle(
            &Vector3::y_axis(),
            f32::consts::FRAC_PI_2,
        ))
        .rotate(&Rotation3::from_axis_angle(
            &Vector3::x_axis(),
            f32::consts::FRAC_PI_2,
        ))
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
        let mut pxy = xy_plane(RectBounds::new(1.0, 1.0));

        let translation_vec = Vector3::new(1.0, -1.0, 3.0);
        let global_point = Point3::from(translation_vec);

        pxy = pxy.translate(&Translation3::from(translation_vec));

        // point at the origin of the plane should be at the
        // translation point
        let converted_global = pxy.get_global_coords(&Point2::origin());
        assert_relative_eq!(global_point, converted_global);
    }

    #[test]
    /// tests coordinate conversion after translation
    fn global_to_local() {
        let t_vec = Vector3::new(1.0, -1.0, 3.0);

        let mut pxy = xy_plane(RectBounds::new(1.0, 1.0)).translate(&Translation3::from(t_vec));

        let global_point = Point3::from(t_vec);

        // point at the origin of the plane should be at the
        // translation point
        let local_point = pxy.get_local_coords(&global_point);
        assert_relative_eq!(local_point, Point2::origin());
    }

    #[test]
    fn xy_plane_coords() {
        let pl_xy = xy_plane(RectBounds::new(1.0, 1.0));

        // test x-axis points
        let globl_x = Point3::new(1.0, 0.0, 0.0);
        let local_x = Point2::new(1.0, 0.0);

        assert_relative_eq!(globl_x, pl_xy.get_global_coords(&local_x));

        // test y-axis points
        let globl_y = Point3::new(0.0, 1.0, 0.0);
        let local_y = Point2::new(0.0, 1.0);

        assert_relative_eq!(globl_y, pl_xy.get_global_coords(&local_y));

        // test x and y
        let globl_xy = Point3::new(-3.0, -3.0, 0.0);
        let local_xy = Point2::new(-3.0, -3.0);

        assert_relative_eq!(globl_xy, pl_xy.get_global_coords(&local_xy));
    }

    #[test]
    fn zx_plane_coords() {
        let pl_zx = zx_plane(RectBounds::new(1.0, 1.0));

        // z-axis becomes local x-axis -> z point should equal local x point
        let globl_z = Point3::new(2.0, 0.0, -1.0);
        let local_z = Point2::new(-1.0, 2.0);

        assert_relative_eq!(globl_z, pl_zx.get_global_coords(&local_z));
    }

    #[test]
    fn yz_plane_coords() {
        let pl_yz = yz_plane(RectBounds::new(1.0, 1.0));

        // y-axis becomes local x-axis -> y point should equal local x point
        let globl_y = Point3::new(0.0, 2.0, -1.0);
        let local_y = Point2::new(2.0, -1.0);

        assert_relative_eq!(globl_y, pl_yz.get_global_coords(&local_y));
    }

    /// Transformation tests
    #[test]
    fn reversible_translation() {
        let mut p1 = xy_plane(RectBounds::new(1.0, 2.0));
        let p2 = xy_plane(RectBounds::new(1.0, 2.0));

        // translate p1 there and back and test equality
        assert_eq!(p1, p2);

        p1 = p1
            .translate(&Translation3::new(1.0, 2.0, 3.0))
            .translate(&Translation3::new(-1.0, -2.0, -3.0));

        assert_eq!(p1, p2);
    }

    #[test]
    fn reversible_rotation() {
        // rotate the xy_plane pi/2 radians around y
        // and back to check reversibility in rotation
        let mut p1 = xy_plane(RectBounds::new(3.0, 3.0))
            .rotate(&Rotation3::from_axis_angle(
                &Vector3::y_axis(),
                f32::consts::FRAC_PI_2,
            ))
            .rotate(&Rotation3::from_axis_angle(
                &Vector3::y_axis(),
                -1.0 * f32::consts::FRAC_PI_2,
            ));

        let p2 = xy_plane(RectBounds::new(3.0, 3.0));
        // note: uses approx equal for floating point numbers!
        assert_eq!(p1, p2);
    }

    // test that global from local and back again gives the same vector
    #[test]
    fn reversible_coord_transform() {
        // test with rotation, translation
        let mut p_yz = yz_plane(RectBounds::new(3.0, 3.0))
            .rotate(&Rotation3::from_axis_angle(
                &Vector3::y_axis(),
                f32::consts::FRAC_PI_2,
            ))
            .translate(&Translation3::new(1.0, 2.0, 3.0));

        let p = Point3::new(-1.0, -2.0, -3.0);

        let local_p = p_yz.get_local_coords(&p);

        assert_relative_eq!(p, p_yz.get_global_coords(&local_p));
    }

    #[test]
    fn inverse_eager() {
        // test with rotation, translation
        let mut p_yz = yz_plane(RectBounds::new(3.0, 3.0))
            .rotate(&Rotation3::from_axis_angle(
                &Vector3::y_axis(),
                f32::consts::FRAC_PI_2,
            ))
            .translate(&Translation3::new(1.0, 2.0, 3.0));

        let p = Point3::new(-1.0, -2.0, -3.0);

        let p_optional = p_yz.get_local_coords(&p);
        let p_eager = p_yz.get_local_coords_eager(&p);

        assert_relative_eq!(p_optional, p_eager);
    }
}
