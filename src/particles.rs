#[derive(Debug, Clone, PartialEq)]
pub struct Particles {
    pub center_x: Vec<f32>,
    pub center_y: Vec<f32>,
    pub radius: Vec<f32>,
    pub velocity_x: Vec<f32>,
    pub velocity_y: Vec<f32>,
    pub mass: Vec<f32>,
    pub count: usize,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            center_x: Vec::with_capacity(capacity),
            center_y: Vec::with_capacity(capacity),
            radius: Vec::with_capacity(capacity),
            velocity_x: Vec::with_capacity(capacity),
            velocity_y: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            count: 0,
        }
    }

    pub fn from_particles(particles: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let count = particles.len();
        let mut result = Self::new(count);

        for (cx, cy, r, vx, vy, m) in particles {
            result.center_x.push(cx);
            result.center_y.push(cy);
            result.radius.push(r);
            result.velocity_x.push(vx);
            result.velocity_y.push(vy);
            result.mass.push(m);
        }
        result.count = count;
        result
    }

    pub fn overlap(&self, i: usize, j: usize) -> (bool, f32) {
        let delta_x = self.center_x[i] - self.center_x[j];
        let delta_y = self.center_y[i] - self.center_y[j];
        let distance_sq = delta_x * delta_x + delta_y * delta_y;
        let radius_sum = self.radius[i] + self.radius[j];
        (distance_sq < radius_sum * radius_sum, distance_sq)
    }

    pub fn clear(&mut self) {
        self.center_x.clear();
        self.center_y.clear();
        self.radius.clear();
        self.velocity_x.clear();
        self.velocity_y.clear();
        self.mass.clear();
        self.count = 0;
    }

    pub fn append(&mut self, other: Self) {
        self.center_x.extend(other.center_x);
        self.center_y.extend(other.center_y);
        self.radius.extend(other.radius);
        self.velocity_x.extend(other.velocity_x);
        self.velocity_y.extend(other.velocity_y);
        self.mass.extend(other.mass);
        self.count += other.count;
    }

    pub fn push(&mut self, particle: (f32, f32, f32, f32, f32, f32)) {
        self.center_x.push(particle.0);
        self.center_y.push(particle.1);
        self.radius.push(particle.2);
        self.velocity_x.push(particle.3);
        self.velocity_y.push(particle.4);
        self.mass.push(particle.5);
        self.count += 1;
    }
}
