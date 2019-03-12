This file has working notes and a "living" implementation plan
that can change as work goes on.

Main goal: describe 3D planes using Rust

Usage notes:
- We choose to not consume transformation objects as they are applied
- However, Bounds objects are owned by their associated Plane 
- Bounds, transformations are all strongly typed

----------------------------------------------------------------

Requirements:
- Data structure for planes X 
- Convert from global coordinates to 2D planar coordinates X
- Data structures for 2D, 3D vectors X (used nalgebra structs)
- Transformation matrices for these vectors handling translation and rotation. X 
  -> Note that a general transform matrix for a 3D vector is a 4x4 w/
     homogenous coordinates.
- All matrices must be invertible X (enforced by nalgebra types)
- Test whether a given point is inside plane bounds X 

Assumptions:
- Planes are all parallel
- Rect bounds are OK for now, but keep genericity in mind
- Assume f32 (c++ double) is OK -> taken from variant_data type

Plan:
- Use nalgebra crate for linear algebra ops, maybe hide this impl detail with
    wrapper classes.
- 3D planes -> use a struct with conversion operators
- Take bounds logic from ACTS impl

Improvements:
- Add tolerance considerations 
- Add unit system (currently unitless)
- Flesh out differences between Detector and Plane Surface (?)
- Add Serde integration for for saving sets of useful planes for different runs

----------------------------------------------------------------

Misc:
- Think about global / local enum for coords?
- This implementation is kind of shaking out like a very
    strongly typed module, making transforms very explicit.
    This has advantages (intent) but is less flexible (e.g.
    have to pass a Translation3 instead of Vector3 to translate.
    This has the effect of exposing the implementation details
    to the caller, but I think the safety and intent advantages
    are enough to outweigh this concern.
