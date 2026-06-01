// Switch the simulation between f32 and f64 by editing the three lines below.
// The GPU rendering path always uses f32 (GLSL `float` / `glVertexAttribPointer`)
// and downcasts at the boundary.

pub type Real = f32;
pub type AtomicReal = atomic_float::AtomicF32;
pub const TAU: Real = std::f64::consts::TAU as Real;
