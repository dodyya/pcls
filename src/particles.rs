use crate::types::{AtomicReal, Real};
use std::sync::atomic::Ordering::Relaxed as O;

#[derive(Debug)]
pub struct Particle {
    x: AtomicReal,
    y: AtomicReal,
    ox: AtomicReal,
    oy: AtomicReal,
    ax: AtomicReal,
    ay: AtomicReal,
    r: AtomicReal,
    m: AtomicReal,
    hue: AtomicReal,
}

impl Particle {
    pub fn new(x: Real, y: Real, r: Real, m: Real, hue: Real) -> Self {
        Self {
            x: AtomicReal::new(x),
            y: AtomicReal::new(y),
            ox: AtomicReal::new(x),
            oy: AtomicReal::new(y),
            ax: AtomicReal::new(0.0),
            ay: AtomicReal::new(0.0),
            r: AtomicReal::new(r),
            m: AtomicReal::new(m),
            hue: AtomicReal::new(hue),
        }
    }
}

#[derive(Debug)]
pub struct Particles {
    pub particles: Vec<Particle>,
    pub count: usize,
    pub g_toward_center: bool,
    pub hue_force_enabled: bool,
    pub donut_enabled: bool,
    // Tracked here so the hue-force can clamp its distance floor to the smallest
    // live particle (mitigates the 1/d² blow-up when particles overlap).
    pub min_radius: Real,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            particles: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
            hue_force_enabled: false,
            donut_enabled: true,
            min_radius: Real::INFINITY,
        }
    }

    pub fn clear(&mut self) {
        self.particles = vec![];
        self.count = 0;
        self.min_radius = Real::INFINITY;
    }

    pub fn push(&mut self, particle: (Real, Real, Real, Real, Real)) -> usize {
        self.particles.push(Particle::new(
            particle.0, particle.1, particle.2, particle.3, particle.4,
        ));
        self.count += 1;
        self.min_radius = self.min_radius.min(particle.2);
        self.count - 1
    }

    pub fn get_x(&self, i: usize) -> Real {
        self.particles[i].x.load(O)
    }
    pub fn get_y(&self, i: usize) -> Real {
        self.particles[i].y.load(O)
    }
    pub fn get_ox(&self, i: usize) -> Real {
        self.particles[i].ox.load(O)
    }
    pub fn get_oy(&self, i: usize) -> Real {
        self.particles[i].oy.load(O)
    }
    pub fn get_ax(&self, i: usize) -> Real {
        self.particles[i].ax.load(O)
    }
    pub fn get_ay(&self, i: usize) -> Real {
        self.particles[i].ay.load(O)
    }
    pub fn get_r(&self, i: usize) -> Real {
        self.particles[i].r.load(O)
    }
    pub fn get_m(&self, i: usize) -> Real {
        self.particles[i].m.load(O)
    }
    pub fn get_hue(&self, i: usize) -> Real {
        self.particles[i].hue.load(O)
    }

    pub fn set_x(&self, i: usize, data: Real) {
        self.particles[i].x.store(data, O)
    }
    pub fn set_y(&self, i: usize, data: Real) {
        self.particles[i].y.store(data, O)
    }
    pub fn set_ox(&self, i: usize, data: Real) {
        self.particles[i].ox.store(data, O)
    }
    pub fn set_oy(&self, i: usize, data: Real) {
        self.particles[i].oy.store(data, O)
    }
    pub fn set_ax(&self, i: usize, data: Real) {
        self.particles[i].ax.store(data, O)
    }
    pub fn set_ay(&self, i: usize, data: Real) {
        self.particles[i].ay.store(data, O)
    }
    // Atomic accumulators so parallel force passes don't lose contributions
    // to the load/store race that plain set_ax(get_ax + d) would have.
    pub fn add_ax(&self, i: usize, delta: Real) {
        self.particles[i].ax.fetch_add(delta, O);
    }
    pub fn add_ay(&self, i: usize, delta: Real) {
        self.particles[i].ay.fetch_add(delta, O);
    }

    pub fn get_drawable(&self) -> impl Iterator<Item = (Real, Real, Real, Real)> + '_ {
        self.particles
            .iter()
            .map(|p| (p.x.load(O), p.y.load(O), p.r.load(O), p.hue.load(O)))
    }
}
