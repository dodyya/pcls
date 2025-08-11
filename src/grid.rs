use crate::util::Array3D;

use crate::particles::{ParticleID, Particles};
const DEPTH: usize = 4;

// pub type GridKey = (i32, i32);

pub struct HashGrid {
    pub cell_count: usize,
    pub map: Array3D<Option<ParticleID>>,
}

impl HashGrid {
    pub fn new(cell_size: f32) -> Self {
        let cell_count = (2.0 / cell_size).ceil() as usize;
        Self {
            cell_count,
            map: Array3D::default(cell_count, cell_count, DEPTH),
        }
    }

    pub fn index(&self, x: f32, y: f32) -> (usize, usize) {
        let cell_size = 2.0 / self.cell_count as f32;
        let x_index = ((x + 1.0) / cell_size).floor() as usize;
        let y_index = ((y + 1.0) / cell_size).floor() as usize;
        (x_index, y_index)
    }

    pub fn try_insert(&mut self, id: usize, x: f32, y: f32) {
        let ind = self.index(x, y);
        let cell = &mut self.map[ind];
        for j in 0..DEPTH {
            if cell[j] == None {
                cell[j] = Some(id);
            }
        }
    }

    pub fn update(&mut self, particles: &Particles) {
        self.map.clear();
        for i in 0..particles.count {
            let ind = self.index(particles.x[i], particles.y[i]);
            let cell = &mut self.map[ind];
            for d in 0..DEPTH {
                if cell[d] == None {
                    cell[d] = Some(i);
                    break;
                }
            }
        }
    }
}
