use crate::grid::{GridKey, HashGrid};
use crate::particles::ParticleID;

use crate::particles::Particles;
pub struct Phx {
    pub pcls: Particles,
    pub grid: HashGrid,
}
const DT: f32 = 1.0 / 12.0;

impl Phx {
    pub fn new(cell_size: f32) -> Self {
        let mut pcls = Particles::new(20_000);
        pcls.add_10k(0.0, 0.0, 0.005, 1.0);
        pcls.add_10k(0.0, 0.0, 0.005, 1.0);
        let mut grid = HashGrid::new(cell_size);
        for i in 0..pcls.count {
            grid.insert(i, pcls.x[i], pcls.y[i]);
        }
        Self { pcls, grid }
    }

    pub fn resolve_overlaps(&mut self) {
        let all_keys: Vec<GridKey> = self.grid.map.keys().cloned().collect();
        let relative_positions = [(1, 0), (1, 1), (0, 1), (-1, 1)];
        for key in all_keys {
            if let Some(my_ids) = self.grid.map.get(&key).cloned() {
                for (idx, &i) in my_ids.iter().enumerate() {
                    for &j in my_ids.iter().skip(idx + 1) {
                        self.pcls.overlap(i, j);
                    }

                    for (dx, dy) in relative_positions.iter() {
                        let cell_key = (key.0 + dx, key.1 + dy);
                        if let Some(indices) = self.grid.map.get(&cell_key) {
                            for &j in indices.iter() {
                                self.pcls.overlap(i, j);
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn step(&mut self) {
        for _ in 0..10 {
            self.pcls.apply_gravity();
            self.pcls.constrain();
            self.grid.update(&self.pcls);
            self.resolve_overlaps();
            self.pcls.verlet(DT / 10.0);
        }
    }

    pub fn get_drawable_particles(&self) -> (&[f32], &[f32], &[f32]) {
        (&self.pcls.x, &self.pcls.y, &self.pcls.radius)
    }

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        self.pcls.push((x, y, radius, vx, vy, mass));
        self.grid.insert(self.pcls.count - 1, x, y);
    }

    pub fn clear(&mut self) {
        self.pcls.clear();
        self.grid.clear();
    }
    pub fn toggle_gravity(&mut self) {
        self.pcls.g_toward_center = !self.pcls.g_toward_center;
    }
}
