use nalgebra::{Affine3, Point3, Rotation3, Translation3, Vector2, Vector3};

/// # A detector geometry module
/// In production, this name would have to be changed to avoid possible
/// name clashes.

/// The Plane data structure, somewhat analagous to PlaneSurface in ACTS
#[derive(Debug)]
pub struct Plane {
    pub centroid: Point3<f32>, // global coords
    pub normal: Vector3<f32>,  // global coords
    pub bounds: (f32, f32),    // plane coords
    pub global_to_local: Affine3<f32> // from the world to the plane
}

impl Plane {
    pub fn get_local_coords(&self, global_vec: &Vector3<f32>) -> Vector2<f32> {
        // broken impl for now
        Vector2::new(0.0, 0.0)
    }

    /// expects translation vector in global coords
    /// -> mutates the current global_to_local transform
    pub fn translate(&mut self, translation_vec: &Translation3<f32>) {
        //
    }

    /// expects rotation matrix in global coords
    /// -> mutates the current global_to_local transform
    pub fn rotate(&mut self, rotation_mat : &Rotation3<f32>) {
        //
    }
}
