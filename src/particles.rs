use atomic_float::AtomicF32;
use std::{
    iter::repeat_n,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc, RwLock,
    },
};

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
const GRAVITY: f32 = 0.1;
// const GRAVITY_TOWARDS_CENTER: bool = false;
const WASHING_MACHINE: bool = false;
const RESTITUTION: f32 = 0.6;
pub(crate) const O: Ordering = Ordering::SeqCst;

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

    pub fn push(&mut self, particle: (f32, f32, f32, f32, f32, f32)) {
        self.x.push(AtomicF32::new(particle.0));
        self.y.push(AtomicF32::new(particle.1));
        self.r.push(AtomicF32::new(particle.2));
        self.ox.push(AtomicF32::new(particle.0));
        self.oy.push(AtomicF32::new(particle.1));
        self.ax.push(AtomicF32::new(0.0));
        self.ay.push(AtomicF32::new(0.0));
        self.m.push(AtomicF32::new(particle.5));

        self.count += 1;
    }

    pub fn apply_gravity(&self) {
        for i in 0..self.count {
            if self.g_toward_center {
                let x = self.x[i].load(O);
                let y = self.y[i].load(O);
                let r2 = x.abs().powi(2) + y.abs().powi(2);
                let v_x = x / (r2.sqrt());
                let v_y = y / (r2.sqrt());
                self.ax[i].store(-GRAVITY * v_x * (1.0 / (r2 + 0.2)), O);
                self.ay[i].store(-GRAVITY * v_y * (1.0 / (r2 + 0.2)), O);
            } else {
                self.ay[i].store(-GRAVITY, O);
            }
        }
    }

    // #[inline(never)]
    pub fn overlap(&self, i: usize, j: usize) {
        let xi = self.x[i].load(O);
        let xj = self.x[j].load(O);
        let yi = self.y[i].load(O);
        let yj = self.y[j].load(O);
        let dx = xi - xj;
        let dy = yi - yj;
        let distance_sq = dx * dx + dy * dy;
        let ri = self.r[i].load(O);
        let rj = self.r[j].load(O);
        if distance_sq > (ri + rj) * (ri + rj) {
            return;
        }

        let distance = distance_sq.sqrt();
        let overlap_distance = ri + rj - distance;
        let (normal_x, normal_y) = if distance != 0.0 {
            (dx / distance, dy / distance)
        } else {
            let theta = rand::thread_rng().gen_range(0.0..std::f32::consts::PI * 2.0);
            (theta.cos(), theta.sin())
        };

        let mi = self.m[i].load(O);
        let mj = self.m[j].load(O);
        let mass_ratio_1 = mi / (mi + mj);
        let mass_ratio_2 = mj / (mi + mj);

        let correction = RESTITUTION * overlap_distance * 0.5;
        let correction_x = correction * normal_x;
        let correction_y = correction * normal_y;

        self.x[i].fetch_add(correction_x * mass_ratio_2, O);
        self.y[i].fetch_add(correction_y * mass_ratio_2, O);
        self.x[j].fetch_sub(correction_x * mass_ratio_1, O);
        self.y[j].fetch_sub(correction_y * mass_ratio_1, O);
    }

    pub fn verlet(&self, dt: f32) {
        for i in 0..self.count {
            let x = self.x[i].load(O);
            let y = self.y[i].load(O);
            let vx = x - self.ox[i].load(O);
            let vy = y - self.oy[i].load(O);
            self.ox[i].store(x, O);
            self.oy[i].store(y, O);
            self.x[i].store(x + vx + self.ax[i].load(O) * dt * dt, O);
            self.y[i].store(y + vy + self.ay[i].load(O) * dt * dt, O);
            self.ax[i].store(0.0, O);
            self.ay[i].store(0.0, O);
        }
    }

    #[inline(never)]
    pub fn constrain(&self) {
        for i in 0..self.count {
            let x = self.x[i].load(O);
            let y = self.y[i].load(O);
            let r = self.r[i].load(O);
            if WASHING_MACHINE {
                let center_dist = (x * x + y * y).sqrt();
                let factor = if center_dist + r > 1.0 {
                    (1.0 - r) / center_dist
                } else if center_dist - r < 0.3 {
                    (0.3 + r) / center_dist
                } else {
                    1.0
                };
                self.x[i].store(x * factor, O);
                self.y[i].store(y * factor, O);
            } else {
                if x + r > 1.0 {
                    self.x[i].store(1.0 - r, O);
                } else if x - r < -1.0 {
                    self.x[i].store(-1.0 + r, O);
                }
                if y + r > 1.0 {
                    self.y[i].store(1.0 - r, O);
                } else if y - r < -1.0 {
                    self.y[i].store(-1.0 + r, O);
                }
            }
        }
    }
    pub fn stop(&self) {
        for i in 0..self.count {
            self.oy[i].store(self.x[i].load(O), O);
            self.oy[i].store(self.y[i].load(O), O);
        }
    }
}
