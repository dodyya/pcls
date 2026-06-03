use crate::types::{AtomicReal, Real};
use std::sync::atomic::Ordering::Relaxed as O;
use std::sync::atomic::AtomicUsize;

// Sentinel for "this particle is not currently registered in any grid cell".
// Stored in cell_xs (cell_ys is irrelevant when unregistered).
pub const CELL_UNREGISTERED: usize = usize::MAX;

// SoA: one Vec per field so sequential reads (e.g. all xs) coalesce in cache
// instead of striding through a 7-atomic struct.
#[derive(Debug)]
pub struct Particles {
    pub xs: Vec<AtomicReal>,
    pub ys: Vec<AtomicReal>,
    pub oxs: Vec<AtomicReal>,
    pub oys: Vec<AtomicReal>,
    pub axs: Vec<AtomicReal>,
    pub ays: Vec<AtomicReal>,
    pub hues: Vec<AtomicReal>,
    // Grid cell each particle is currently registered in. CELL_UNREGISTERED
    // means "not in the grid right now" (e.g. fresh particle pre-insert, or
    // last insert hit a full cell).
    pub cell_xs: Vec<AtomicUsize>,
    pub cell_ys: Vec<AtomicUsize>,
    pub count: usize,
    pub g_toward_center: bool,
    pub hue_force_enabled: bool,
    pub donut_enabled: bool,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            xs: Vec::with_capacity(capacity),
            ys: Vec::with_capacity(capacity),
            oxs: Vec::with_capacity(capacity),
            oys: Vec::with_capacity(capacity),
            axs: Vec::with_capacity(capacity),
            ays: Vec::with_capacity(capacity),
            hues: Vec::with_capacity(capacity),
            cell_xs: Vec::with_capacity(capacity),
            cell_ys: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
            hue_force_enabled: false,
            donut_enabled: true,
        }
    }

    pub fn clear(&mut self) {
        self.xs.clear();
        self.ys.clear();
        self.oxs.clear();
        self.oys.clear();
        self.axs.clear();
        self.ays.clear();
        self.hues.clear();
        self.cell_xs.clear();
        self.cell_ys.clear();
        self.count = 0;
    }

    pub fn push(&mut self, particle: (Real, Real, Real)) -> usize {
        let (x, y, hue) = particle;
        self.xs.push(AtomicReal::new(x));
        self.ys.push(AtomicReal::new(y));
        self.oxs.push(AtomicReal::new(x));
        self.oys.push(AtomicReal::new(y));
        self.axs.push(AtomicReal::new(0.0));
        self.ays.push(AtomicReal::new(0.0));
        self.hues.push(AtomicReal::new(hue));
        self.cell_xs.push(AtomicUsize::new(CELL_UNREGISTERED));
        self.cell_ys.push(AtomicUsize::new(0));
        self.count += 1;
        self.count - 1
    }

    pub fn get_x(&self, i: usize) -> Real {
        self.xs[i].load(O)
    }
    pub fn get_y(&self, i: usize) -> Real {
        self.ys[i].load(O)
    }
    pub fn get_ox(&self, i: usize) -> Real {
        self.oxs[i].load(O)
    }
    pub fn get_oy(&self, i: usize) -> Real {
        self.oys[i].load(O)
    }
    pub fn get_ax(&self, i: usize) -> Real {
        self.axs[i].load(O)
    }
    pub fn get_ay(&self, i: usize) -> Real {
        self.ays[i].load(O)
    }
    pub fn get_hue(&self, i: usize) -> Real {
        self.hues[i].load(O)
    }

    pub fn set_x(&self, i: usize, data: Real) {
        self.xs[i].store(data, O)
    }
    pub fn set_y(&self, i: usize, data: Real) {
        self.ys[i].store(data, O)
    }
    pub fn set_ox(&self, i: usize, data: Real) {
        self.oxs[i].store(data, O)
    }
    pub fn set_oy(&self, i: usize, data: Real) {
        self.oys[i].store(data, O)
    }
    pub fn set_ax(&self, i: usize, data: Real) {
        self.axs[i].store(data, O)
    }
    pub fn set_ay(&self, i: usize, data: Real) {
        self.ays[i].store(data, O)
    }
    // Atomic accumulators so parallel force passes don't lose contributions
    // to the load/store race that plain set_ax(get_ax + d) would have.
    pub fn add_ax(&self, i: usize, delta: Real) {
        self.axs[i].fetch_add(delta, O);
    }
    pub fn add_ay(&self, i: usize, delta: Real) {
        self.ays[i].fetch_add(delta, O);
    }

    // Currently-registered grid cell, or CELL_UNREGISTERED if not in the grid.
    pub fn get_cell(&self, i: usize) -> (usize, usize) {
        (self.cell_xs[i].load(O), self.cell_ys[i].load(O))
    }
    pub fn set_cell(&self, i: usize, cell: (usize, usize)) {
        self.cell_xs[i].store(cell.0, O);
        self.cell_ys[i].store(cell.1, O);
    }
    pub fn set_cell_unregistered(&self, i: usize) {
        self.cell_xs[i].store(CELL_UNREGISTERED, O);
    }

    pub fn get_drawable(&self) -> impl Iterator<Item = (Real, Real, Real)> + '_ {
        (0..self.count).map(move |i| (self.xs[i].load(O), self.ys[i].load(O), self.hues[i].load(O)))
    }
}
