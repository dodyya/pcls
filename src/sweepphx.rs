use crate::particles::{ParticleID, Particles};
use crate::sweep::SweepAndPrune;

pub struct Phx {
    pub pcls: Particles,
    pub broadphase: SweepAndPrune,
}

const DT: f32 = 1.0 / 12.0;
impl Phx {
    pub fn new(cell_size: f32) -> Self {
        let mut particles = Particles::new(20_000);
        particles.add_10k(0.0, 0.0, 0.005, 1.0);
        particles.add_10k(0.0, 0.0, 0.005, 1.0);
        let broadphase = SweepAndPrune::new(&particles);

        Self {
            pcls: particles,
            broadphase,
        }
    }

    pub fn resolve_overlaps(&mut self) {
        let entries = &self.broadphase.entries_x;

        for i in 0..entries.len() {
            let entry_i = entries[i];

            for j in (i + 1)..entries.len() {
                let entry_j = entries[j];

                if entry_j.min > entry_i.max {
                    break;
                }
                self.pcls.overlap(entry_i.id, entry_j.id);
            }
        }
        self.broadphase.update(&self.pcls);
    }

    pub fn step(&mut self) {
        for _ in 0..10 {
            self.pcls.apply_gravity();
            self.pcls.constrain();
            self.broadphase.update(&self.pcls);
            self.resolve_overlaps();
            self.pcls.verlet(DT / 10.0);
        }
    }

    pub fn get_drawable_particles(&self) -> (&[f32], &[f32], &[f32]) {
        (&self.pcls.x, &self.pcls.y, &self.pcls.radius)
    }

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        self.pcls.push((x, y, radius, vx, vy, mass));
        self.broadphase.insert(self.pcls.count - 1, x, y, radius);
    }

    pub fn clear(&mut self) {
        self.pcls.clear();
        self.broadphase.clear();
    }

    pub fn toggle_gravity(&mut self) {
        self.pcls.g_toward_center = !self.pcls.g_toward_center;
    }
}
