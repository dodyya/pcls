use crate::types::Real;
use cuda_core::device_buffer::DeviceCopy;

// Cold device-side integration state (16 B). Split out from positions so that
// neighbour-gather kernels (hue_force, overlap, grid_build) fetch only the 12 B
// of XyHue they actually need instead of dragging this along.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct RawParticle {
    pub ox: Real,
    pub oy: Real,
    pub ax: Real,
    pub ay: Real,
}

// Safety: RawParticle is repr(C) over Copy primitives only.
unsafe impl DeviceCopy for RawParticle {}

// Hot position record (12 B): the (x, y, hue) read by hue_force, overlap, grid_build,
// gravity, constrain, verlet, and the renderer. The `extract` step copies this buffer
// straight into the CUDA-mapped VBO (same layout, same stride), and the GL VAO binds
// to it as (a_pos vec2 @0, a_hue float @8, stride 12).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct XyHue {
    pub x: Real,
    pub y: Real,
    pub hue: Real,
}

// Safety: XyHue is repr(C) over Copy primitives only.
unsafe impl DeviceCopy for XyHue {}

// Host-side staging mirror: positions + integration state held as parallel vectors.
// Written by add_particle, read on htod, never updated from device except during a
// grow refresh. The renderer never touches it.
pub struct Particles {
    pub positions: Vec<XyHue>,
    pub integs: Vec<RawParticle>,
    pub count: usize,
    pub g_toward_center: bool,
    pub hue_force_enabled: bool,
    pub donut_enabled: bool,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            positions: Vec::with_capacity(capacity),
            integs: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
            hue_force_enabled: false,
            donut_enabled: true,
        }
    }

    pub fn clear(&mut self) {
        self.positions.clear();
        self.integs.clear();
        self.count = 0;
    }

    pub fn push(&mut self, x: Real, y: Real, hue: Real) {
        self.positions.push(XyHue { x, y, hue });
        self.integs.push(RawParticle {
            ox: x,
            oy: y,
            ax: 0.0,
            ay: 0.0,
        });
        self.count += 1;
    }
}
