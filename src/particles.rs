use crate::types::{AtomicReal, Real};
use std::sync::atomic::Ordering::Relaxed as O;
use std::sync::atomic::AtomicUsize;

// Sentinel for "this particle is not currently registered in any grid cell".
// Stored in Particle::cell_x (cell_y is irrelevant when unregistered).
pub const CELL_UNREGISTERED: usize = usize::MAX;

// AoS layout (one Vec<Particle>, all fields per particle contiguous). Beat the
// SoA layout on bench: donut, gravity-down, hue on, 12 substeps, 60 measured
// steps after 40 warmup, best-of-3, total ms/step:
//   N      SoA    AoS     Δ
//   1000   3.66   3.66    ~0%
//   3000   4.51   4.42   −2.0%
//   6000   5.63   5.48   −2.7%
//   10000  7.33   7.01   −4.2%
// Per-phase: verlet dropped 17–29% (the win — it touches 6 atomics per particle
// and AoS fits one particle in ~one cache line). Hue/overlap moved within noise
// (pair-based loads see the same total bytes either way).
#[derive(Debug, Default)]
pub struct Particle {
    pub x: AtomicReal,
    pub y: AtomicReal,
    pub ox: AtomicReal,
    pub oy: AtomicReal,
    pub ax: AtomicReal,
    pub ay: AtomicReal,
    pub hue: AtomicReal,
    // Grid cell each particle is currently registered in. CELL_UNREGISTERED in
    // cell_x means "not in the grid right now" (fresh, or last insert hit a
    // full cell).
    pub cell_x: AtomicUsize,
    pub cell_y: AtomicUsize,
}

#[derive(Debug)]
pub struct Particles {
    pub items: Vec<Particle>,
    pub count: usize,
    pub g_toward_center: bool,
    pub hue_force_enabled: bool,
    pub donut_enabled: bool,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            items: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
            hue_force_enabled: false,
            donut_enabled: true,
        }
    }

    pub fn clear(&mut self) {
        self.items.clear();
        self.count = 0;
    }

    pub fn push(&mut self, particle: (Real, Real, Real)) -> usize {
        let (x, y, hue) = particle;
        self.items.push(Particle {
            x: AtomicReal::new(x),
            y: AtomicReal::new(y),
            ox: AtomicReal::new(x),
            oy: AtomicReal::new(y),
            ax: AtomicReal::new(0.0),
            ay: AtomicReal::new(0.0),
            hue: AtomicReal::new(hue),
            cell_x: AtomicUsize::new(CELL_UNREGISTERED),
            cell_y: AtomicUsize::new(0),
        });
        self.count += 1;
        self.count - 1
    }

    pub fn get_x(&self, i: usize) -> Real { self.items[i].x.load(O) }
    pub fn get_y(&self, i: usize) -> Real { self.items[i].y.load(O) }
    pub fn get_ox(&self, i: usize) -> Real { self.items[i].ox.load(O) }
    pub fn get_oy(&self, i: usize) -> Real { self.items[i].oy.load(O) }
    pub fn get_ax(&self, i: usize) -> Real { self.items[i].ax.load(O) }
    pub fn get_ay(&self, i: usize) -> Real { self.items[i].ay.load(O) }
    pub fn get_hue(&self, i: usize) -> Real { self.items[i].hue.load(O) }

    pub fn set_x(&self, i: usize, data: Real) { self.items[i].x.store(data, O) }
    pub fn set_y(&self, i: usize, data: Real) { self.items[i].y.store(data, O) }
    pub fn set_ox(&self, i: usize, data: Real) { self.items[i].ox.store(data, O) }
    pub fn set_oy(&self, i: usize, data: Real) { self.items[i].oy.store(data, O) }
    pub fn set_ax(&self, i: usize, data: Real) { self.items[i].ax.store(data, O) }
    pub fn set_ay(&self, i: usize, data: Real) { self.items[i].ay.store(data, O) }

    // Atomic accumulators so parallel force passes don't lose contributions
    // to the load/store race that plain set_ax(get_ax + d) would have.
    pub fn add_ax(&self, i: usize, delta: Real) { self.items[i].ax.fetch_add(delta, O); }
    pub fn add_ay(&self, i: usize, delta: Real) { self.items[i].ay.fetch_add(delta, O); }

    pub fn get_cell(&self, i: usize) -> (usize, usize) {
        (self.items[i].cell_x.load(O), self.items[i].cell_y.load(O))
    }
    pub fn set_cell(&self, i: usize, cell: (usize, usize)) {
        self.items[i].cell_x.store(cell.0, O);
        self.items[i].cell_y.store(cell.1, O);
    }
    pub fn set_cell_unregistered(&self, i: usize) {
        self.items[i].cell_x.store(CELL_UNREGISTERED, O);
    }

    pub fn get_drawable(&self) -> impl Iterator<Item = (Real, Real, Real)> + '_ {
        (0..self.count).map(move |i| {
            let p = &self.items[i];
            (p.x.load(O), p.y.load(O), p.hue.load(O))
        })
    }
}
