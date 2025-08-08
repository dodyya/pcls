#[derive(Debug, Clone, PartialEq)]
pub struct Particles {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub radius: Vec<f32>,
    pub vx: Vec<f32>,
    pub vy: Vec<f32>,
    pub mass: Vec<f32>,
    pub count: usize,
}
const RESTITUTION: f32 = 0.5;
const VELOCITY_DAMPING: f32 = 1.00;
const COLLISION_DAMPING: f32 = 1.00;
const GRAVITY: f32 = 0.01;
const WALL_DAMPING: f32 = 0.9;
const GRAVITY_TOWARDS_CENTER: bool = false;

pub type ParticleID = usize;

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            radius: Vec::with_capacity(capacity),
            vx: Vec::with_capacity(capacity),
            vy: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            count: 0,
        }
    }

    pub fn from_particles(particles: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let count = particles.len();
        let mut result = Self::new(count);

        for (cx, cy, r, vx, vy, m) in particles {
            result.x.push(cx);
            result.y.push(cy);
            result.radius.push(r);
            result.vx.push(vx);
            result.vy.push(vy);
            result.mass.push(m);
        }
        result.count = count;
        result
    }

    pub fn overlap_info(&self, i: usize, j: usize) -> (bool, f32) {
        let delta_x = self.x[i] - self.x[j];
        let delta_y = self.y[i] - self.y[j];
        let distance_sq = delta_x * delta_x + delta_y * delta_y;
        let radius_sum = self.radius[i] + self.radius[j];
        (distance_sq < radius_sum * radius_sum, distance_sq)
    }

    pub fn clear(&mut self) {
        self.x.clear();
        self.y.clear();
        self.radius.clear();
        self.vx.clear();
        self.vy.clear();
        self.mass.clear();
        self.count = 0;
    }

    pub fn push(&mut self, particle: (f32, f32, f32, f32, f32, f32)) {
        self.x.push(particle.0);
        self.y.push(particle.1);
        self.radius.push(particle.2);
        self.vx.push(particle.3);
        self.vy.push(particle.4);
        self.mass.push(particle.5);
        self.count += 1;
    }

    pub fn collision(&mut self, i: usize, j: usize) {
        let delta_x = self.x[i] - self.x[j];
        let delta_y = self.y[i] - self.y[j];
        let distance_sq = delta_x * delta_x + delta_y * delta_y;
        let distance = distance_sq.sqrt();

        let (normal_x, normal_y) = if distance != 0.0 {
            (delta_x / distance, delta_y / distance)
        } else {
            (1.0, 0.0)
        };

        let relative_velocity_x = self.vx[i] - self.vx[j];
        let relative_velocity_y = self.vy[i] - self.vy[j];
        let velocity_along_normal = relative_velocity_x * normal_x + relative_velocity_y * normal_y;

        if velocity_along_normal > 0.0 {
            return;
        }

        let a_inv = 1.0 / self.mass[i];
        let b_inv = 1.0 / self.mass[j];
        let impulse_mag = -(1.0 + RESTITUTION) * velocity_along_normal / (a_inv + b_inv);
        let impulse_x = normal_x * impulse_mag * COLLISION_DAMPING;
        let impulse_y = normal_y * impulse_mag * COLLISION_DAMPING;

        self.vx[i] += impulse_x * a_inv;
        self.vy[i] += impulse_y * a_inv;
        self.vx[j] -= impulse_x * b_inv;
        self.vy[j] -= impulse_y * b_inv;
    }
    /// Apply velocity damping to all particles.
    pub fn apply_velocity_damping(&mut self) {
        for i in 0..self.count {
            self.vx[i] *= VELOCITY_DAMPING;
            self.vy[i] *= VELOCITY_DAMPING;
        }
    }
    /// Apply gravity to all particles.
    pub fn apply_gravity(&mut self) {
        for i in 0..self.count {
            if GRAVITY_TOWARDS_CENTER {
                self.vx[i] -= GRAVITY * self.x[i];
                self.vy[i] -= GRAVITY * self.y[i];
            } else {
                self.vy[i] -= 0.2 * GRAVITY;
            }
        }
    }

    pub fn overlap(&mut self, i: usize, j: usize) -> bool {
        if i == j {
            return false;
        }

        let (overlap, distance_sq) = self.overlap_info(i, j);
        if !overlap {
            return false;
        }

        let distance = distance_sq.sqrt();
        let delta_x = self.x[i] - self.x[j];
        let delta_y = self.y[i] - self.y[j];

        let (normal_x, normal_y) = if distance != 0.0 {
            (delta_x / distance, delta_y / distance)
        } else {
            (1.0, 0.0)
        };

        let overlap_distance = self.radius[i] + self.radius[j] - distance;
        let correction_x = normal_x * (overlap_distance * 0.5);
        let correction_y = normal_y * (overlap_distance * 0.5);

        let old_pos_i_x = self.x[i];
        let old_pos_i_y = self.y[i];
        let old_pos_j_x = self.x[j];
        let old_pos_j_y = self.y[j];

        self.x[i] += correction_x;
        self.y[i] += correction_y;
        self.x[j] -= correction_x;
        self.y[j] -= correction_y;
        true
    }

    pub fn update_kinematics(&mut self) {
        for i in 0..self.count {
            self.x[i] += self.vx[i];
            self.y[i] += self.vy[i];
        }
    }

    pub fn resolve_wall_collisions(&mut self) {
        for i in 0..self.count {
            let radius = self.radius[i];

            if self.x[i] - radius < -1.0 {
                self.x[i] = -1.0 + radius;
                self.vx[i] = WALL_DAMPING * self.vx[i].abs();
            }

            if self.x[i] + radius > 1.0 {
                self.x[i] = 1.0 - radius;
                self.vx[i] = -WALL_DAMPING * self.vx[i].abs();
            }

            if self.y[i] - radius < -1.0 {
                self.y[i] = -1.0 + radius;
                self.vy[i] = WALL_DAMPING * self.vy[i].abs();
            }

            if self.y[i] + radius > 1.0 {
                self.y[i] = 1.0 - radius;
                self.vy[i] = -WALL_DAMPING * self.vy[i].abs();
            }
        }
    }
}
