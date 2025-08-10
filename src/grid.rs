use crate::util::Array2D;
use std::collections::HashMap;

use crate::particles::{ParticleID, Particles};

// pub type GridKey = (i32, i32);

pub struct HashGrid {
    pub cell_count: usize,
    // pub map: HashMap<GridKey, Vec<ParticleID>>,
    pub map: Array2D<Vec<ParticleID>>,
}

impl HashGrid {
    pub fn new(cell_size: f32) -> Self {
        let cell_count = (2.0 / cell_size).ceil() as usize;
        Self {
            cell_count,
            map: Array2D::default(cell_count, cell_count),
        }
    }

    // pub fn key(&self, x: f32, y: f32) -> GridKey {
    //     (
    //         (x / 2.0 * self.cell_count as f32).floor() as i32,
    //         (y / 2.0 * self.cell_count as f32).floor() as i32,
    //     )
    // }

    pub fn index(&self, x: f32, y: f32) -> (usize, usize) {
        let cell_size = 2.0 / self.cell_count as f32;
        let x_index = ((x + 1.0) / cell_size).floor() as usize;
        let y_index = ((y + 1.0) / cell_size).floor() as usize;
        (x_index, y_index)
    }

    pub fn insert(&mut self, id: usize, x: f32, y: f32) {
        let i = self.index(x, y);
        self.map[i].push(id);
    }

    pub fn update(&mut self, particles: &Particles) {
        self.clear();
        for i in 0..particles.count {
            let ind = self.index(particles.x[i], particles.y[i]);
            self.map[ind].push(i);
        }
    }

    pub fn clear(&mut self) {
        self.map
            .data
            .iter_mut()
            .for_each(|cell_vec| cell_vec.clear());
    }
}
