This file has working notes and a "living" implementation plan
that can change as work goes on.

Main goal: describe 3D planes using Rust

Usage notes:
- We choose to not consume transformation objects as they are applied
- Bounds, transformations are all strongly typed

----------------------------------------------------------------

Requirements:
- Data structure for planes
- Convert from global coordinates to 2D planar coordinates
- Data structures for 2D, 3D vectors
- Transformation matrices for these vectors handling translation and rotation.
  -> Note that a general transform matrix for a 3D vector is a 4x4 w/
     homogenous coordinates.
- All matrices must be invertible
- Test whether a given point is inside plane bounds

Assumptions:
- Planes are all parallel
- Plane bounds are OK for now, but keep genericity in mind
- Assume f32 (c++ double) is OK -> taken from variant_data type

Plan:
- Use nalgebra crate for linear algebra ops, maybe hide this impl detail with
    wrapper classes.
- 3D planes -> use a struct with conversion operators
- only calculate inversions when they are asked for since matrix inversion is relatively expensive (use Option type)
- Take bounds logic from ACTS impl
- if we have time, use serde for saving sets of useful planes for different runs

----------------------------------------------------------------

Misc:
- Think about global / local enum for coords?
- This implementation is kind of shaking out like a very
    strongly typed module, making transforms very explicit.
    This has advantages (intent) but is less flexible (e.g.
    have to pass a Translation3 instead of Vector3 to translate.
    This has the effect of exposing the implementation details
    to the caller, but I think the safety and intent advantages
    are enough to outweight this concern.