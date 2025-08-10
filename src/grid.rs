use std::collections::HashMap;

use crate::particles::{ParticleID, Particles};

pub type GridKey = (i32, i32);

pub struct HashGrid {
    cell_count: i32,
    pub map: HashMap<GridKey, Vec<ParticleID>>,
}

impl HashGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_count: (2.0 / cell_size).ceil() as i32,
            map: HashMap::new(),
        }
    }

    pub fn key(&self, x: f32, y: f32) -> GridKey {
        (
            (x / 2.0 * self.cell_count as f32).floor() as i32,
            (y / 2.0 * self.cell_count as f32).floor() as i32,
        )
    }

    pub fn insert(&mut self, id: usize, x: f32, y: f32) {
        let grid_key = self.key(x, y);
        self.map.entry(grid_key).or_default().push(id);
    }

    pub fn update(&mut self, particles: &Particles) {
        self.map.clear();
        for i in 0..particles.count {
            self.map
                .entry(self.key(particles.x[i], particles.y[i]))
                .or_default()
                .push(i);
        }
    }

    pub fn clear(&mut self) {
        self.map.clear();
    }
}
