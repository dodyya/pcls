use atomic_float::AtomicF32;
use std::sync::atomic::Ordering::Relaxed as O;

#[derive(Debug)]
pub struct Particle {
    x: AtomicF32,
    y: AtomicF32,
    ox: AtomicF32,
    oy: AtomicF32,
    ax: AtomicF32,
    ay: AtomicF32,
    r: AtomicF32,
    m: AtomicF32,
    charge: AtomicF32,
}

impl Particle {
    pub fn new(x: f32, y: f32, r: f32, m: f32, charge: f32) -> Self {
        Self {
            x: AtomicF32::new(x),
            y: AtomicF32::new(y),
            ox: AtomicF32::new(x),
            oy: AtomicF32::new(y),
            ax: AtomicF32::new(0.0),
            ay: AtomicF32::new(0.0),
            r: AtomicF32::new(r),
            m: AtomicF32::new(m),
            charge: AtomicF32::new(charge),
        }
    }
}

#[derive(Debug)]
pub struct Particles {
    pub particles: Vec<Particle>,
    pub count: usize,
    pub g_toward_center: bool,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            particles: Vec::with_capacity(capacity),
            count: 0,
            g_toward_center: false,
        }
    }

    pub fn add_10k(&mut self, x: f32, y: f32, r: f32, m: f32, charge: f32) {
        self.particles.reserve(10_000);
        for _ in 0..10_000 {
            self.particles.push(Particle::new(x, y, r, m, charge));
        }
        self.count += 10_000;
    }

    pub fn clear(&mut self) {
        self.particles = vec![];
        self.count = 0;
    }

    pub fn push(&mut self, particle: (f32, f32, f32, f32, f32)) -> usize {
        self.particles.push(Particle::new(
            particle.0, particle.1, particle.2, particle.3, particle.4,
        ));
        self.count += 1;
        self.count - 1
    }

    pub fn get_x(&self, i: usize) -> f32 {
        self.particles[i].x.load(O)
    }
    pub fn get_y(&self, i: usize) -> f32 {
        self.particles[i].y.load(O)
    }
    pub fn get_ox(&self, i: usize) -> f32 {
        self.particles[i].ox.load(O)
    }
    pub fn get_oy(&self, i: usize) -> f32 {
        self.particles[i].oy.load(O)
    }
    pub fn get_ax(&self, i: usize) -> f32 {
        self.particles[i].ax.load(O)
    }
    pub fn get_ay(&self, i: usize) -> f32 {
        self.particles[i].ay.load(O)
    }
    pub fn get_r(&self, i: usize) -> f32 {
        self.particles[i].r.load(O)
    }
    pub fn get_m(&self, i: usize) -> f32 {
        self.particles[i].m.load(O)
    }
    pub fn get_c(&self, i: usize) -> f32 {
        self.particles[i].charge.load(O)
    }

    pub fn set_x(&self, i: usize, data: f32) {
        self.particles[i].x.store(data, O)
    }
    pub fn set_y(&self, i: usize, data: f32) {
        self.particles[i].y.store(data, O)
    }
    pub fn set_ox(&self, i: usize, data: f32) {
        self.particles[i].ox.store(data, O)
    }
    pub fn set_oy(&self, i: usize, data: f32) {
        self.particles[i].oy.store(data, O)
    }
    pub fn set_ax(&self, i: usize, data: f32) {
        self.particles[i].ax.store(data, O)
    }
    pub fn set_ay(&self, i: usize, data: f32) {
        self.particles[i].ay.store(data, O)
    }

    pub fn get_drawable(&self) -> impl Iterator<Item = (f32, f32, f32, f32)> + '_ {
        self.particles
            .iter()
            .map(|p| (p.x.load(O), p.y.load(O), p.r.load(O), p.charge.load(O)))
    }
}
