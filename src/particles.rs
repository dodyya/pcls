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
}
const GRAVITY: f32 = 0.1;
const GRAVITY_TOWARDS_CENTER: bool = true;
const WASHING_MACHINE: bool = false;
const RESTITUTION: f32 = 0.75;

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
        }
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

    pub fn apply_gravity(&mut self) {
        for i in 0..self.count {
            if GRAVITY_TOWARDS_CENTER {
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

    pub fn overlap(&mut self, i: usize, j: usize) -> bool {
        if i == j {
            return false;
        }
        let dx = self.x[i] - self.x[j];
        let dy = self.y[i] - self.y[j];
        let distance_sq = dx * dx + dy * dy;
        let distance = distance_sq.sqrt();
        if distance > (self.radius[i] + self.radius[j]) {
            return false;
        }

        let (normal_x, normal_y) = if distance != 0.0 {
            (dx / distance, dy / distance)
        } else {
            (1.0, 0.0)
        };

        let mass_ratio_1 = self.mass[i] / (self.mass[i] + self.mass[j]);
        let mass_ratio_2 = self.mass[j] / (self.mass[i] + self.mass[j]);

        let overlap_distance = self.radius[i] + self.radius[j] - distance;
        let correction_x = RESTITUTION * normal_x * overlap_distance * 0.5;
        let correction_y = RESTITUTION * normal_y * overlap_distance * 0.5;

        self.x[i] += correction_x * mass_ratio_2;
        self.y[i] += correction_y * mass_ratio_2;
        self.x[j] -= correction_x * mass_ratio_1;
        self.y[j] -= correction_y * mass_ratio_1;
        true
    }

    pub fn verlet(&mut self, dt: f32) {
        for i in 0..self.count {
            let vx = self.x[i] - self.ox[i];
            let vy = self.y[i] - self.oy[i];
            self.ox[i] = self.x[i];
            self.oy[i] = self.y[i];
            self.x[i] = self.x[i] + vx + self.ax[i] * dt * dt;
            self.y[i] = self.y[i] + vy + self.ay[i] * dt * dt;
            self.ax[i] = 0.0;
            self.ay[i] = 0.0;
        }
    }

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
                if self.x[i] + self.radius[i] > 1.0 {
                    self.x[i] = 1.0 - self.radius[i];
                } else if self.x[i] - self.radius[i] < -1.0 {
                    self.x[i] = -1.0 + self.radius[i];
                }
                if self.y[i] + self.radius[i] > 1.0 {
                    self.y[i] = 1.0 - self.radius[i];
                } else if self.y[i] - self.radius[i] < -1.0 {
                    self.y[i] = -1.0 + self.radius[i];
                }
            }
        }
    }
}
