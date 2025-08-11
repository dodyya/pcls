use rand::Rng;
#[derive(Debug, Clone, PartialEq)]
pub struct Particles {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub ox: Vec<f32>,
    pub oy: Vec<f32>,
    pub ax: Vec<f32>,
    pub ay: Vec<f32>,
    pub radius: Vec<f32>,
    pub mass: Vec<f32>,
    pub count: usize,
    pub g_toward_center: bool,
}
const GRAVITY: f32 = 0.1;
// const GRAVITY_TOWARDS_CENTER: bool = false;
const WASHING_MACHINE: bool = false;
const RESTITUTION: f32 = 0.6;

pub type ParticleID = usize;

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            ox: Vec::with_capacity(capacity),
            oy: Vec::with_capacity(capacity),
            ax: Vec::with_capacity(capacity),
            ay: Vec::with_capacity(capacity),
            radius: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
        }
    }

    pub fn add_10k(&mut self, x: f32, y: f32, r: f32, m: f32) {
        for _ in 0..10_000 {
            self.x.push(x);
            self.y.push(y);
            self.radius.push(r);
            self.ox.push(x);
            self.oy.push(y);
            self.ax.push(0.0);
            self.ay.push(0.0);
            self.mass.push(m);
        }
        self.count += 10_000;
    }

    pub fn from_particles(particles: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let count = particles.len();
        let mut result = Self::new(count);

        for (cx, cy, r, ox, oy, m) in particles {
            result.x.push(cx);
            result.y.push(cy);
            result.radius.push(r);
            result.ox.push(ox);
            result.oy.push(oy);
            result.mass.push(m);
        }
        result.count = count;
        result
    }

    pub fn clear(&mut self) {
        self.x.clear();
        self.y.clear();
        self.radius.clear();
        self.ox.clear();
        self.oy.clear();
        self.mass.clear();
        self.count = 0;
    }

    pub fn push(&mut self, particle: (f32, f32, f32, f32, f32, f32)) {
        self.x.push(particle.0);
        self.y.push(particle.1);
        self.radius.push(particle.2);
        self.ox.push(particle.0);
        self.oy.push(particle.1);
        self.ax.push(0.0);
        self.ay.push(0.0);
        self.mass.push(particle.5);
        self.count += 1;
    }

    #[inline(never)]
    pub fn apply_gravity(&mut self) {
        for i in 0..self.count {
            if self.g_toward_center {
                let r2 = self.x[i].abs().powi(2) + self.y[i].abs().powi(2);
                let v_x = self.x[i] / (r2.sqrt());
                let v_y = self.y[i] / (r2.sqrt());
                self.ax[i] = -GRAVITY * v_x * (1.0 / (r2 + 0.2));
                self.ay[i] = -GRAVITY * v_y * (1.0 / (r2 + 0.2));
            } else {
                self.ay[i] = -GRAVITY;
            }
        }
    }

    // #[inline(never)]
    pub unsafe fn overlap(&mut self, i: usize, j: usize) {
        let &xi = self.x.get_unchecked(i);
        let &xj = self.x.get_unchecked(j);
        let &yi = self.y.get_unchecked(i);
        let &yj = self.y.get_unchecked(j);
        let dx = xi - xj;
        let dy = yi - yj;
        let distance_sq = dx * dx + dy * dy;
        let ri = self.radius.get_unchecked(i);
        let rj = self.radius.get_unchecked(j);
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

        let mi = self.mass.get_unchecked(i);
        let mj = self.mass.get_unchecked(j);
        let mass_ratio_1 = mi / (mi + mj);
        let mass_ratio_2 = mj / (mi + mj);

        let correction = RESTITUTION * overlap_distance * 0.5;
        let correction_x = correction * normal_x;
        let correction_y = correction * normal_y;

        *self.x.get_unchecked_mut(i) = xi + correction_x * mass_ratio_2;
        *self.y.get_unchecked_mut(i) = yi + correction_y * mass_ratio_2;
        *self.x.get_unchecked_mut(j) = xj - correction_x * mass_ratio_1;
        *self.y.get_unchecked_mut(j) = yj - correction_y * mass_ratio_1;
    }

    pub fn verlet(&mut self, dt: f32) {
        for i in 0..self.count {
            let x = self.x[i];
            let y = self.y[i];
            let vx = x - self.ox[i];
            let vy = y - self.oy[i];
            self.ox[i] = x;
            self.oy[i] = y;
            self.x[i] = x + vx + self.ax[i] * dt * dt;
            self.y[i] = y + vy + self.ay[i] * dt * dt;
            self.ax[i] = 0.0;
            self.ay[i] = 0.0;
        }
    }

    #[inline(never)]
    pub fn constrain(&mut self) {
        for i in 0..self.count {
            if WASHING_MACHINE {
                let center_dist = (self.x[i] * self.x[i] + self.y[i] * self.y[i]).sqrt();
                if center_dist + self.radius[i] > 1.0 {
                    let factor = (1.0 - self.radius[i]) / center_dist;
                    self.x[i] *= factor;
                    self.y[i] *= factor;
                } else if center_dist - self.radius[i] < 0.3 {
                    let factor = (0.3 + self.radius[i]) / center_dist;
                    self.x[i] *= factor;
                    self.y[i] *= factor;
                }
            } else {
                let r = self.radius[i];
                if self.x[i] + r > 1.0 {
                    self.x[i] = 1.0 - r;
                } else if self.x[i] - r < -1.0 {
                    self.x[i] = -1.0 + r;
                }
                if self.y[i] + r > 1.0 {
                    self.y[i] = 1.0 - r;
                } else if self.y[i] - r < -1.0 {
                    self.y[i] = -1.0 + r;
                }
            }
        }
    }
    pub fn stop(&mut self) {
        self.ox = self.x.clone();
        self.oy = self.y.clone();
    }
}
