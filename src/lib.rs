//#[cfg(test)]
//mod tests {
//    // import all names from the outer scope
//    use super::*;
//
//    /// Transformation tests
//    #[test]
//    fn reversible_translation() {
//
//    }
//
//    #[test]
//    fn reversible_rotation() {
//
//    }
//
//    // test that global from local and back again gives the same vector
//    #[test]
//    fn reversible_coord_transform() {
//        let p1 = Plane {
//            centroid: Point3::new(0.0, 0.0, 0.0),
//            normal: Vector3::new(0.0, 0.0, 1.0),
//            bounds: (1.0, 2.0),
//            global_to_local: Affine3::identity()
//        };
//
//        let p2 = p1.clone();
//
//    }
//
//}
