use crate::particles::{ParticleID, Particles};
#[derive(Debug, Clone, Copy)]
struct AxisEntry {
    id: ParticleID,
    min: f32,
    max: f32,
}

pub struct SweepAndPrune {
    entries_x: Vec<AxisEntry>,
    entries_y: Vec<AxisEntry>,
}

impl SweepAndPrune {
    pub fn new(particles: &Particles) -> Self {
        let mut entries_x = Vec::with_capacity(particles.count);
        let mut entries_y = Vec::with_capacity(particles.count);

        for i in 0..particles.count {
            let radius = particles.radius[i];

            entries_x.push(AxisEntry {
                id: i,
                min: particles.x[i] - radius,
                max: particles.x[i] + radius,
            });

            entries_y.push(AxisEntry {
                id: i,
                min: particles.y[i] - radius,
                max: particles.y[i] + radius,
            });
        }

        entries_x.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        entries_y.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());

        Self {
            entries_x,
            entries_y,
        }
    }

    pub fn clear(&mut self) {
        self.entries_x.clear();
        self.entries_y.clear();
    }

    pub fn update(&mut self, particles: &Particles) {
        for entry in &mut self.entries_x {
            let i = entry.id;
            entry.min = particles.x[i] - particles.radius[i];
            entry.max = particles.x[i] + particles.radius[i];
        }

        for entry in &mut self.entries_y {
            let i = entry.id;
            entry.min = particles.y[i] - particles.radius[i];
            entry.max = particles.y[i] + particles.radius[i];
        }

        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        self.entries_y
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }

    pub fn get_potential_collisions(&self) -> Vec<(usize, usize)> {
        let mut potential_pairs = Vec::new();
        let entries = &self.entries_x;

        for i in 0..entries.len() {
            let entry_i = &entries[i];

            for j in (i + 1)..entries.len() {
                let entry_j = &entries[j];

                if entry_j.min > entry_i.max {
                    break;
                }

                potential_pairs.push((entry_i.id, entry_j.id));
            }
        }

        potential_pairs
    }

    pub fn rebuild(&mut self, particles: &Particles) {
        self.entries_x.clear();
        self.entries_y.clear();

        for i in 0..particles.count {
            let radius = particles.radius[i];

            self.entries_x.push(AxisEntry {
                id: i,
                min: particles.x[i] - radius,
                max: particles.x[i] + radius,
            });

            self.entries_y.push(AxisEntry {
                id: i,
                min: particles.y[i] - radius,
                max: particles.y[i] + radius,
            });
        }

        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        self.entries_y
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }
}
