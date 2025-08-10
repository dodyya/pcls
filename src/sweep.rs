use crate::particles::{ParticleID, Particles};
#[derive(Debug, Clone, Copy)]
pub struct AxisEntry {
    pub id: ParticleID,
    pub min: f32,
    pub max: f32,
}

pub struct SweepAndPrune {
    pub entries_x: Vec<AxisEntry>,
    // pub entries_y: Vec<AxisEntry>,
}

impl SweepAndPrune {
    pub fn new(particles: &Particles) -> Self {
        let mut entries_x = Vec::with_capacity(particles.count);
        // let mut entries_y = Vec::with_capacity(particles.count);

        for i in 0..particles.count {
            let radius = particles.radius[i];

            entries_x.push(AxisEntry {
                id: i,
                min: particles.x[i] - radius,
                max: particles.x[i] + radius,
            });

            // entries_y.push(AxisEntry {
            //     id: i,
            //     min: particles.y[i] - radius,
            //     max: particles.y[i] + radius,
            // });
        }

        entries_x.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        // entries_y.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());

        Self {
            entries_x,
            // entries_y,
        }
    }

    pub fn clear(&mut self) {
        self.entries_x.clear();
        // self.entries_y.clear();
    }

    pub fn update(&mut self, particles: &Particles) {
        for entry in &mut self.entries_x {
            let x = particles.x[entry.id];
            let r = particles.radius[entry.id];
            entry.min = x - r;
            entry.max = x + r;
        }

        // for entry in &mut self.entries_y {
        //     let y = particles.y[entry.id];
        //     let r = particles.radius[entry.id];
        //     entry.min = y - r;
        //     entry.max = y + r;
        // }

        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        // self.entries_y
        //     .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }

    pub fn rebuild(&mut self, particles: &Particles) {
        self.entries_x.clear();
        // self.entries_y.clear();

        for i in 0..particles.count {
            let radius = particles.radius[i];

            self.entries_x.push(AxisEntry {
                id: i,
                min: particles.x[i] - radius,
                max: particles.x[i] + radius,
            });

            // self.entries_y.push(AxisEntry {
            //     id: i,
            //     min: particles.y[i] - radius,
            //     max: particles.y[i] + radius,
            // });
        }

        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        // self.entries_y
        //     .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }

    pub fn insert(&mut self, id: usize, x: f32, y: f32, radius: f32) {
        self.entries_x.push(AxisEntry {
            id,
            min: x - radius,
            max: x + radius,
        });

        // self.entries_y.push(AxisEntry {
        //     id,
        //     min: y - radius,
        //     max: y + radius,
        // });

        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        // self.entries_y
        //     .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }
}
