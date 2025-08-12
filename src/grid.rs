use crate::array::Array3D;

use crate::maybe_id::MaybeID;
use crate::particles::Particles;

#[derive(Debug)]
pub struct Grid {
    pub cell_count: usize,
    pub depth: usize,
    pub map: Array3D<MaybeID>,
}

impl Grid {
    pub fn new(cell_size: f32, depth: usize) -> Self {
        let cell_count = (2.0 / cell_size).ceil() as usize;
        Self {
            cell_count,
            depth,
            map: Array3D::default(cell_count, cell_count, depth),
        }
    }

    pub fn index(&self, x: f32, y: f32) -> (usize, usize) {
        let cell_size = 2.0 / self.cell_count as f32;
        let x_index = ((x + 1.0) / cell_size).floor() as usize;
        let y_index = ((y + 1.0) / cell_size).floor() as usize;
        (x_index, y_index)
    }

    pub fn try_insert(&self, id: usize, x: f32, y: f32) {
        let ind = self.index(x, y);
        let cell = &self.map[ind];
        for j in 0..self.depth {
            if cell[j].is_none() {
                cell[j].set(id);
                break;
            }
        }
    }

    pub fn update(&self, particles: &Particles) {
        self.map.clear();
        for i in 0..particles.count {
            let x = particles.get_x(i);
            let y = particles.get_y(i);
            let ind = self.index(x, y);
            let cell = &self.map[ind];
            for d in 0..self.depth {
                if cell[d].is_none() {
                    cell[d].set(i);
                    break;
                }
            }
        }
    }
}
