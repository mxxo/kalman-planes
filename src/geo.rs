use nalgebra::{Affine3, Point3, Rotation3, Translation3, Vector2, Vector3};

/// # A detector geometry module
/// In production, this name would have to be changed to avoid possible
/// name clashes.

/// The bounds struct is used to distinguish Bounds types from other tuples
#[derive(Debug, PartialEq)]
pub struct Bounds{
    x: f32,
    y: f32,
}

impl Bounds {
    pub fn new(x: f32, y: f32) -> Bounds {
        // rough safety check for zero or negative bounds
        assert!(x > 0.0 && y > 0.0);
        Bounds { x, y }
    }
}

/// The Plane data structure, somewhat analagous to PlaneSurface in ACTS
#[derive(Debug)]
pub struct Plane {
    pub centroid: Point3<f32>, // global coords
    pub normal: Vector3<f32>,  // global coords
    pub bounds: Bounds,    // plane coords
    pub global_to_local: Affine3<f32> // from the world to the plane
}

/// We only compare the geometry and transformation matrix
/// of the two planes i.e. the history of how the plane
/// got to its state doesn't matter

impl PartialEq for Plane {
    fn eq(&self, other: &Plane) -> bool {
        self.centroid == other.centroid
        && self.normal == other.normal
        && self.bounds == other.bounds
        && self.global_to_local == other.global_to_local
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
    pub fn translate(&mut self, translation_vec: Translation3<f32>) {
        self.centroid = translation_vec * self.centroid;
        self.global_to_local = translation_vec * self.global_to_local;
    }

    /// expects rotation matrix in global coords
    /// -> mutates the current global_to_local transform
    pub fn rotate(&mut self, rotation_mat : &Rotation3<f32>) {
        //
    }
}

/// # Convenience constructors for planes
/// TODO need to think about how the transform matrix should be constructed

/// xy plane has a positive z-direction normal
pub fn xy_plane(bounds: Bounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 0.0, 1.0),
        bounds,
        global_to_local: Affine3::identity()
    }
}

/// xz plane has a positive y-direction normal
pub fn xz_plane(bounds: Bounds) -> Plane {
    Plane {
        centroid: Point3::new(0.0, 0.0, 0.0),
        normal: Vector3::new(0.0, 1.0, 0.0),
        bounds,
        global_to_local: Affine3::identity()
    }
}

/// yz plane has a positive x-direction normal
pub fn yz_plane(bounds: Bounds) -> Plane {
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

    /// Bounds tests
    #[test]
    #[should_panic]
    fn positive_bounds() {
       Bounds::new(-1.0, 0.0);
    }

    #[test]
    #[should_panic]
    fn zero_bounds() {
        Bounds::new(0.0, 1.0);
    }

    /// Transformation tests
    #[test]
    fn reversible_translation() {

        let mut p1 = xy_plane( Bounds::new(1.0, 2.0) );
        let p2 = xy_plane( Bounds::new(1.0, 2.0) );

        // translate p1 there and back and test equality
        assert_eq!(p1, p2);

        p1.translate(Translation3::new(1.0, 2.0, 3.0));
        p1.translate(Translation3::new(-1.0, -2.0, -3.0));

        assert_eq!(p1, p2);
    }

    #[test]
    fn reversible_rotation() {
        assert_eq!(1, 2);
    }

    // test that global from local and back again gives the same vector
    #[test]
    fn reversible_coord_transform() {
        //let mut p1 = Plane {
        //    centroid: Point3::new(0.0, 0.0, 0.0),
        //    normal: Vector3::new(0.0, 0.0, 1.0),
        //    bounds: Bounds{ x: 1.0, y: 2.0 },
        //    global_to_local: Affine3::identity()
        //};

        //let p2 = Plane {
        //    centroid: Point3::origin(),
        //    normal: Vector3::new(0.0, 0.0, 1.0),
        //    bounds: Bounds{ x: 1.0, y: 2.0 },
        //    global_to_local: Affine3::identity()
        //};
        assert_eq!(1, 2);
    }
}
