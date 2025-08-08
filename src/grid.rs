use std::collections::HashMap;

use crate::particles::{ParticleID, Particles};

pub type GridKey = (i32, i32);

pub struct HashGrid {
    cell_size: f32,
    cell_count: i32,
    pub map: HashMap<GridKey, Vec<ParticleID>>,
}

impl HashGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cell_count: (2.0 / cell_size).ceil() as i32,
            map: HashMap::new(),
        }
    }

    pub fn get_grid_key(&self, x: f32, y: f32) -> GridKey {
        (
            (x / 2.0 * self.cell_count as f32).floor() as i32,
            (y / 2.0 * self.cell_count as f32).floor() as i32,
        )
    }

    pub fn insert(&mut self, index: usize, x: f32, y: f32) {
        let grid_key = self.get_grid_key(x, y);
        self.map.entry(grid_key).or_default().push(index);
    }

    pub fn remove(&mut self, index: usize, x: f32, y: f32) {
        let grid_key = self.get_grid_key(x, y);
        if let Some(vec) = self.map.get_mut(&grid_key) {
            if let Some(pos) = vec.iter().position(|&i| i == index) {
                vec.swap_remove(pos);
                if vec.is_empty() {
                    self.map.remove(&grid_key);
                }
            }
        }
    }

    pub fn get_possible_neighbors(&self, key: GridKey) -> Vec<usize> {
        let mut neighbors = vec![];

        if let Some(indices) = self.map.get(&key) {
            neighbors.extend(indices.iter().cloned());
        }

        let relative_positions = [(1, 0), (1, 1), (0, 1), (-1, 1)];

        for (dx, dy) in relative_positions.iter() {
            let cell_key = (key.0 + dx, key.1 + dy);
            if let Some(indices) = self.map.get(&cell_key) {
                neighbors.extend(indices.iter());
            }
        }

        neighbors
    }

    pub fn update_particle(
        &mut self,
        index: usize,
        old_x: f32,
        old_y: f32,
        new_x: f32,
        new_y: f32,
    ) {
        let old_key = self.get_grid_key(old_x, old_y);
        let new_key = self.get_grid_key(new_x, new_y);
        if old_key != new_key {
            if let Some(cell) = self.map.get_mut(&old_key) {
                if let Some(pos) = cell.iter().position(|&i| i == index) {
                    cell.swap_remove(pos);
                    if cell.is_empty() {
                        self.map.remove(&old_key);
                    }
                }
            }
            self.map.entry(new_key).or_default().push(index);
        }
    }

    //Just reinitialize
    pub fn update(&mut self, particles: &Particles) {
        self.map.clear();
        for i in 0..particles.count {
            let key = self.get_grid_key(particles.x[i], particles.y[i]);
            self.map.entry(key).or_default().push(i);
        }
    }

    pub fn clear(&mut self) {
        self.map.clear();
    }
}
