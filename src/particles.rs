use atomic_float::AtomicF32;
use std::sync::atomic::Ordering::Relaxed as O;

use rand::Rng;
#[derive(Debug)]
pub struct Particles {
    pub x: Vec<AtomicF32>,
    pub y: Vec<AtomicF32>,
    pub ox: Vec<AtomicF32>,
    pub oy: Vec<AtomicF32>,
    pub ax: Vec<AtomicF32>,
    pub ay: Vec<AtomicF32>,
    pub r: Vec<AtomicF32>,
    pub m: Vec<AtomicF32>,
    pub count: usize,
    pub g_toward_center: bool,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            ox: Vec::with_capacity(capacity),
            oy: Vec::with_capacity(capacity),
            ax: Vec::with_capacity(capacity),
            ay: Vec::with_capacity(capacity),
            r: Vec::with_capacity(capacity),
            m: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
        }
    }

    pub fn add_10k(&mut self, x: f32, y: f32, r: f32, m: f32) {
        for _ in 0..10_000 {
            self.x.push(AtomicF32::new(x));
            self.y.push(AtomicF32::new(y));
            self.r.push(AtomicF32::new(r));
            self.ox.push(AtomicF32::new(x));
            self.oy.push(AtomicF32::new(y));
            self.ax.push(AtomicF32::new(0.0));
            self.ay.push(AtomicF32::new(0.0));
            self.m.push(AtomicF32::new(m));
        }
        self.count += 10_000;
    }

    pub fn clear(&mut self) {
        self.x = vec![];
        self.y = vec![];
        self.r = vec![];
        self.ox = vec![];
        self.oy = vec![];
        self.m = vec![];
        self.count = 0;
    }

    pub fn push(&mut self, particle: (f32, f32, f32, f32, f32, f32)) -> usize {
        self.x.push(AtomicF32::new(particle.0));
        self.y.push(AtomicF32::new(particle.1));
        self.r.push(AtomicF32::new(particle.2));
        self.ox.push(AtomicF32::new(particle.0));
        self.oy.push(AtomicF32::new(particle.1));
        self.ax.push(AtomicF32::new(0.0));
        self.ay.push(AtomicF32::new(0.0));
        self.m.push(AtomicF32::new(particle.5));

        self.count += 1;
        self.count - 1
    }

    pub fn get_x(&self, i: usize) -> f32 {
        self.x[i].load(O)
    }
    pub fn get_y(&self, i: usize) -> f32 {
        self.y[i].load(O)
    }
    pub fn get_ox(&self, i: usize) -> f32 {
        self.ox[i].load(O)
    }
    pub fn get_oy(&self, i: usize) -> f32 {
        self.oy[i].load(O)
    }
    pub fn get_ax(&self, i: usize) -> f32 {
        self.ax[i].load(O)
    }
    pub fn get_ay(&self, i: usize) -> f32 {
        self.ay[i].load(O)
    }
    pub fn get_r(&self, i: usize) -> f32 {
        self.r[i].load(O)
    }
    pub fn get_m(&self, i: usize) -> f32 {
        self.m[i].load(O)
    }

    pub fn set_x(&self, i: usize, data: f32) {
        self.x[i].store(data, O)
    }
    pub fn set_y(&self, i: usize, data: f32) {
        self.y[i].store(data, O)
    }
    pub fn set_ox(&self, i: usize, data: f32) {
        self.ox[i].store(data, O)
    }
    pub fn set_oy(&self, i: usize, data: f32) {
        self.oy[i].store(data, O)
    }
    pub fn set_ax(&self, i: usize, data: f32) {
        self.ax[i].store(data, O)
    }
    pub fn set_ay(&self, i: usize, data: f32) {
        self.ay[i].store(data, O)
    }
}
