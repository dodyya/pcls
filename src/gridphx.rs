use crate::grid::{GridKey, HashGrid};

use crate::particles::Particles;
pub struct Phx {
    pub pcls: Particles,
    pub grid: HashGrid,
}

const WALL_DAMPING: f32 = 0.9;

impl Phx {
    pub fn new(cell_size: f32, p_data: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let pcls = Particles::from_particles(p_data);
        let mut grid = HashGrid::new(cell_size);
        for i in 0..pcls.count {
            grid.insert(i, pcls.x[i], pcls.y[i]);
        }
        Self { pcls, grid }
    }

    // pub fn update_kinematics(&mut self) {
    //     for i in 0..self.pcls.count {
    //         let old_x = self.pcls.x[i];
    //         let old_y = self.pcls.y[i];

    //         self.pcls.x[i] += self.pcls.vx[i];
    //         self.pcls.y[i] += self.pcls.vy[i];

    //         self.grid
    //             .update_particle(i, old_x, old_y, self.pcls.x[i], self.pcls.y[i]);
    //     }
    // }

    // pub fn resolve_wall_collisions(&mut self) {
    //     for i in 0..self.pcls.count {
    //         let old_x = self.pcls.x[i];
    //         let old_y = self.pcls.y[i];
    //         let radius = self.pcls.radius[i];

    //         if self.pcls.x[i] - radius < -1.0 {
    //             self.pcls.x[i] = -1.0 + radius;
    //             self.pcls.vx[i] = WALL_DAMPING * self.pcls.vx[i].abs();
    //         }
    //         if self.pcls.x[i] + radius > 1.0 {
    //             self.pcls.x[i] = 1.0 - radius;
    //             self.pcls.vx[i] = -WALL_DAMPING * self.pcls.vx[i].abs();
    //         }
    //         if self.pcls.y[i] - radius < -1.0 {
    //             self.pcls.y[i] = -1.0 + radius;
    //             self.pcls.vy[i] = WALL_DAMPING * self.pcls.vy[i].abs();
    //         }
    //         if self.pcls.y[i] + radius > 1.0 {
    //             self.pcls.y[i] = 1.0 - radius;
    //             self.pcls.vy[i] = -WALL_DAMPING * self.pcls.vy[i].abs();
    //         }

    //         self.grid
    //             .update_particle(i, old_x, old_y, self.pcls.x[i], self.pcls.y[i]);
    //     }
    // }

    pub fn resolve_collisions(&mut self) {
        let keys: Vec<GridKey> = self.grid.map.keys().cloned().collect();

        for cell_key in keys {
            let all_relevant_indices = self.grid.get_possible_neighbors(cell_key);

            if let Some(current_cell_indices) = self.grid.map.get(&cell_key).cloned() {
                for (idx, &i) in current_cell_indices.iter().enumerate() {
                    for &j in current_cell_indices.iter().skip(idx + 1) {
                        let (overlap, _) = self.pcls.overlap_info(i, j);
                        if overlap {
                            self.pcls.collision(i, j);
                        }
                    }

                    for &j in &all_relevant_indices {
                        if current_cell_indices.contains(&j) {
                            continue;
                        }

                        let (overlap, _) = self.pcls.overlap_info(i, j);
                        if overlap {
                            self.pcls.collision(i, j);
                        }
                    }
                }
            }
        }
    }

    fn resolve_overlap_between(&mut self, i: usize, j: usize) -> bool {
        if !self.pcls.overlap(i, j) {
            return false;
        }
        let (old_pos_i_x, old_pos_i_y) = (self.pcls.x[i], self.pcls.y[i]);
        let (old_pos_j_x, old_pos_j_y) = (self.pcls.x[j], self.pcls.y[j]);

        self.grid
            .update_particle(i, old_pos_i_x, old_pos_i_y, self.pcls.x[i], self.pcls.y[i]);
        self.grid
            .update_particle(j, old_pos_j_x, old_pos_j_y, self.pcls.x[j], self.pcls.y[j]);
        true
    }

    /// Iteratively resolve overlaps until no collisions are detected or until max_iterations is reached.
    pub fn resolve_overlaps(&mut self, max_iterations: usize) {
        for _ in 0..max_iterations {
            let mut collision_found = false;

            let cell_keys: Vec<GridKey> = self.grid.map.keys().cloned().collect();

            for key in cell_keys {
                let all_relevant_indices = self.grid.get_possible_neighbors(key);

                if let Some(current_cell_indices) = self.grid.map.get(&key).cloned() {
                    for (idx, &i) in current_cell_indices.iter().enumerate() {
                        for &j in current_cell_indices.iter().skip(idx + 1) {
                            if self.resolve_overlap_between(i, j) {
                                collision_found = true;
                            }
                        }

                        for &j in &all_relevant_indices {
                            if current_cell_indices.contains(&j) {
                                continue;
                            }

                            if self.resolve_overlap_between(i, j) {
                                collision_found = true;
                            }
                        }
                    }
                }
            }

            if !collision_found {
                break;
            }
        }
    }

    /// Update the simulation state.
    pub fn step(&mut self) {
        self.pcls.update_kinematics();
        self.pcls.resolve_wall_collisions();
        self.resolve_collisions();
        self.pcls.apply_velocity_damping();
        self.pcls.apply_gravity();
        self.resolve_overlaps(3);
    }

    pub fn get_drawable_particles(&self) -> Vec<(f32, f32, f32)> {
        let mut result = Vec::with_capacity(self.pcls.count);
        for i in 0..self.pcls.count {
            result.push((self.pcls.x[i], self.pcls.y[i], self.pcls.radius[i]));
        }
        result
    }

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        self.pcls.push((x, y, radius, vx, vy, mass));
        self.grid.insert(self.pcls.count - 1, x, y);
    }

    /// Clear all particles and grid data
    pub fn clear(&mut self) {
        self.pcls.clear();
        self.grid.clear();
    }
}
