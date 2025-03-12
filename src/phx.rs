use dashmap::DashMap;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::{
    collections::HashMap,
    io::stdin,
    ops::{Add, AddAssign, Div, Mul, Sub, SubAssign},
    sync::Mutex,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Particles {
    pub center_x: Vec<f32>,
    pub center_y: Vec<f32>,
    pub radius: Vec<f32>,
    pub velocity_x: Vec<f32>,
    pub velocity_y: Vec<f32>,
    pub mass: Vec<f32>,
}

impl Particles {
    fn overlap(&self, i1: usize, i2: usize) -> (bool, f32) {
        let dx = self.center_x[i1] - self.center_x[i2];
        let dy = self.center_y[i1] - self.center_y[i2];
        let distance_sq = dx * dx + dy * dy;
        let radius_sum = self.radius[i1] + self.radius[i2];
        (distance_sq < radius_sum * radius_sum, distance_sq)
    }

    pub fn clear(&mut self) {
        self.center_x.clear();
        self.center_y.clear();
        self.radius.clear();
        self.velocity_x.clear();
        self.velocity_y.clear();
        self.mass.clear();
    }

    pub fn push(
        &mut self,
        center_x: f32,
        center_y: f32,
        radius: f32,
        velocity_x: f32,
        velocity_y: f32,
        mass: f32,
    ) {
        self.center_x.push(center_x);
        self.center_y.push(center_y);
        self.radius.push(radius);
        self.velocity_x.push(velocity_x);
        self.velocity_y.push(velocity_y);
        self.mass.push(mass);
    }

    pub fn len(&self) -> usize {
        self.center_x.len()
    }
}
// Constant wall damping coefficient:
const WALL_DAMPING: f32 = 0.9;
const VELOCITY_DAMPING: f32 = 1.00;
const RESTITUTION: f32 = 0.5;
const GRAVITY: f32 = -0.001;
const COLLISION_DAMPING: f32 = 1.00;

pub type GridKey = (i32, i32);

/// The grid now stores indices (into the Phx.particles Vec)
pub struct HashGrid {
    cell_size: f32,
    cell_count: i32,
    grid: DashMap<GridKey, Vec<usize>>,
}

impl HashGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cell_count: (2.0 / cell_size).ceil() as i32,
            grid: DashMap::new(),
        }
    }

    pub fn get_grid_key(&self, x: f32, y: f32) -> GridKey {
        (
            (x / 2.0 * self.cell_count as f32).floor() as i32,
            (y / 2.0 * self.cell_count as f32).floor() as i32,
        )
    }

    /// Inserts a particle index into the cell determined by p’s position.
    pub fn insert(&mut self, index: usize, x: f32, y: f32) {
        let grid_key = self.get_grid_key(x, y);
        self.grid.entry(grid_key).or_default().push(index);
    }

    /// Removes a given particle index from the grid cell corresponding to a given particle position.
    pub fn remove(&mut self, index: usize, x: f32, y: f32) {
        let grid_key = self.get_grid_key(x, y);
        if let Some(mut vec) = self.grid.get_mut(&grid_key) {
            if let Some(pos) = vec.iter().position(|&i| i == index) {
                vec.swap_remove(pos);
            }
        }
    }

    /// Returns a vector of indices from the current cell and all neighbors.
    pub fn get_possible_neighbors(&self, key: GridKey) -> Vec<usize> {
        let mut neighbors = vec![];
        for i in -1..=1 {
            for j in -1..=1 {
                let cell_key = (key.0 + i, key.1 + j);
                if let Some(indices) = self.grid.get(&cell_key) {
                    neighbors.extend(indices.iter().cloned());
                }
            }
        }
        neighbors
    }

    pub fn clear(&mut self) {
        self.grid.clear();
    }
}

/// Central simulation structure: a particle store and a hash grid.
pub struct Phx {
    pub particles: Particles,
    pub grid: HashGrid,
}

// Add these before impl Phx:
unsafe impl Sync for Phx {}
unsafe impl Send for Phx {}

impl Phx {
    /// Create a new simulation with a given grid cell size.
    pub fn new(cell_size: f32) -> Self {
        let grid = HashGrid::new(cell_size);
        let particles = Particles {
            center_x: vec![],
            center_y: vec![],
            radius: vec![],
            velocity_x: vec![],
            velocity_y: vec![],
            mass: vec![],
        };

        // Insert all particle indices into the grid.
        Self { particles, grid }
    }

    pub fn push(
        &mut self,
        center_x: f32,
        center_y: f32,
        radius: f32,
        velocity_x: f32,
        velocity_y: f32,
        mass: f32,
    ) {
        self.particles
            .push(center_x, center_y, radius, velocity_x, velocity_y, mass);
        let new_index = self.particles.len() - 1;
        self.grid.insert(new_index, center_x, center_y);
    }

    /// Update kinematics: update positions and grid membership.
    pub fn update_kinematics(&mut self) {
        // Clone old positions for grid membership check
        let old_x = self.particles.center_x.clone();
        let old_y = self.particles.center_y.clone();

        // Update positions in parallel
        self.particles
            .center_x
            .par_iter_mut()
            .zip(self.particles.velocity_x.par_iter())
            .for_each(|(p, v)| *p += v);

        self.particles
            .center_y
            .par_iter_mut()
            .zip(self.particles.velocity_y.par_iter())
            .for_each(|(p, v)| *p += v);

        self.update_grid_membership(old_x, old_y);
    }

    fn update_grid_membership(&self, old_x: Vec<f32>, old_y: Vec<f32>) {
        let updates: Vec<(usize, GridKey, GridKey)> = self
            .particles
            .center_x
            .par_iter()
            .zip(self.particles.center_y.par_iter())
            .enumerate()
            .map(|(i, (&x, &y))| {
                let old_key = self.grid.get_grid_key(old_x[i], old_y[i]);
                let new_key = self.grid.get_grid_key(x, y);
                (i, old_key, new_key)
            })
            .filter(|(_, old_key, new_key)| old_key != new_key)
            .collect();

        updates.par_iter().for_each(|&(i, old_key, new_key)| {
            if let Some(mut cell) = self.grid.grid.get_mut(&old_key) {
                cell.retain(|&idx| idx != i);
                if cell.is_empty() {
                    drop(cell);
                    self.grid.grid.remove(&old_key);
                }
            }
            self.grid.grid.entry(new_key).or_default().push(i);
        })
    }

    /// Resolve collisions with walls.
    pub fn resolve_wall_collisions(&mut self) {
        let old_x = self.particles.center_x.clone();
        let old_y = self.particles.center_y.clone();

        for i in 0..self.particles.center_x.len() {
            let r = self.particles.radius[i];
            if self.particles.center_x[i] - r < -1.0 {
                self.particles.center_x[i] = -1.0 + r;
                self.particles.velocity_x[i] = WALL_DAMPING * self.particles.velocity_x[i].abs();
            }
            if self.particles.center_x[i] + r > 1.0 {
                self.particles.center_x[i] = 1.0 - r;
                self.particles.velocity_x[i] = -WALL_DAMPING * self.particles.velocity_x[i].abs();
            }
            if self.particles.center_y[i] - r < -1.0 {
                self.particles.center_y[i] = -1.0 + r;
                self.particles.velocity_y[i] = WALL_DAMPING * self.particles.velocity_y[i].abs();
            }
            if self.particles.center_y[i] + r > 1.0 {
                self.particles.center_y[i] = 1.0 - r;
                self.particles.velocity_y[i] = -WALL_DAMPING * self.particles.velocity_y[i].abs();
            }
        }

        self.update_grid_membership(old_x, old_y);
    }

    /// Resolve collisions between particles.
    pub fn resolve_collisions(&mut self) {
        // Get grid keys to iterate over.
        let keys: Vec<GridKey> = self.grid.grid.iter().map(|entry| *entry.key()).collect();

        let mut collision_pairs = Vec::new();
        for cell_key in keys {
            let neighbor_indices = self.grid.get_possible_neighbors(cell_key);
            // For each particle in the current cell, check collisions with all neighbors.
            if let Some(cell_indices_ref) = self.grid.grid.get(&cell_key) {
                let cell_indices: Vec<usize> = cell_indices_ref.value().clone();
                for &i in &cell_indices {
                    // For each neighbor index (avoid self‑collision)
                    for &j in &neighbor_indices {
                        if i <= j {
                            continue;
                        }
                        // Resolve collision between particles[i] and particles[j]
                        let (overlap, _) = self.particles.overlap(i, j);
                        if overlap {
                            collision_pairs.push((i, j));
                        }
                    }
                }
            }
        }

        for (i, j) in collision_pairs {
            self.resolve_collision_between(i, j);
        }
    }

    fn resolve_collision_between(&mut self, i: usize, j: usize) {
        let dx = self.particles.center_x[i] - self.particles.center_x[j];
        let dy = self.particles.center_y[i] - self.particles.center_y[j];
        let distance_sq = dx * dx + dy * dy;
        let distance = distance_sq.sqrt();

        let normal_x = if distance != 0.0 { dx / distance } else { 1.0 };
        let normal_y = if distance != 0.0 { dy / distance } else { 0.0 };

        // Relative velocity along the normal.
        let relative_vx = self.particles.velocity_x[i] - self.particles.velocity_x[j];
        let relative_vy = self.particles.velocity_y[i] - self.particles.velocity_y[j];
        let velocity_along_normal = relative_vx * normal_x + relative_vy * normal_y;
        if velocity_along_normal > 0.0 {
            return;
        }

        let a_inv = 1.0 / self.particles.mass[i];
        let b_inv = 1.0 / self.particles.mass[j];
        let impulse_mag = -(1.0 + RESTITUTION) * velocity_along_normal / (a_inv + b_inv);
        let imp_x = normal_x * impulse_mag * COLLISION_DAMPING;
        let imp_y = normal_y * impulse_mag * COLLISION_DAMPING;

        self.particles.velocity_x[i] += imp_x * a_inv;
        self.particles.velocity_y[i] += imp_y * a_inv;
        self.particles.velocity_x[j] -= imp_x * b_inv;
        self.particles.velocity_y[j] -= imp_y * b_inv;
    }

    /// Apply velocity damping to all particles.
    pub fn apply_velocity_damping(&mut self) {
        self.particles
            .velocity_x
            .par_iter_mut()
            .for_each(|v| *v *= VELOCITY_DAMPING);
        self.particles
            .velocity_y
            .par_iter_mut()
            .for_each(|v| *v *= VELOCITY_DAMPING);
    }

    fn resolve_overlap_between(&mut self, i: usize, j: usize) -> bool {
        // Avoid self-collision.
        if i == j {
            return false;
        }

        // Compute the squared distance between particle centers.
        let dx = self.particles.center_x[i] - self.particles.center_x[j];
        let dy = self.particles.center_y[i] - self.particles.center_y[j];
        let distance_sq = dx * dx + dy * dy;
        let radius_sum = self.particles.radius[i] + self.particles.radius[j];

        // If there’s no overlap, exit early.
        if distance_sq >= radius_sum * radius_sum {
            return false;
        }

        let distance = distance_sq.sqrt();
        let (normal_x, normal_y) = if distance != 0.0 {
            (dx / distance, dy / distance)
        } else {
            (1.0, 0.0)
        };

        // Determine how much to correct each particle.
        let overlap_distance = radius_sum - distance;
        let correction_x = normal_x * (overlap_distance * 0.5);
        let correction_y = normal_y * (overlap_distance * 0.5);

        // Save old positions for grid update.
        let old_x_i = self.particles.center_x[i];
        let old_y_i = self.particles.center_y[i];
        let old_x_j = self.particles.center_x[j];
        let old_y_j = self.particles.center_y[j];

        // Adjust positions.
        self.particles.center_x[i] += correction_x;
        self.particles.center_y[i] += correction_y;
        self.particles.center_x[j] -= correction_x;
        self.particles.center_y[j] -= correction_y;

        // Update grid membership for particle i.
        let old_key_i = self.grid.get_grid_key(old_x_i, old_y_i);
        let new_key_i = self
            .grid
            .get_grid_key(self.particles.center_x[i], self.particles.center_y[i]);
        if old_key_i != new_key_i {
            if let Some(mut cell) = self.grid.grid.get_mut(&old_key_i) {
                cell.retain(|&idx| idx != i);
                if cell.is_empty() {
                    drop(cell);
                    self.grid.grid.remove(&old_key_i);
                }
            }
            self.grid.grid.entry(new_key_i).or_default().push(i);
        }

        // Update grid membership for particle j.
        let old_key_j = self.grid.get_grid_key(old_x_j, old_y_j);
        let new_key_j = self
            .grid
            .get_grid_key(self.particles.center_x[j], self.particles.center_y[j]);
        if old_key_j != new_key_j {
            if let Some(mut cell) = self.grid.grid.get_mut(&old_key_j) {
                cell.retain(|&idx| idx != j);
                if cell.is_empty() {
                    drop(cell);
                    self.grid.grid.remove(&old_key_j);
                }
            }
            self.grid.grid.entry(new_key_j).or_default().push(j);
        }

        true
    }

    /// Resolve overlaps using a 9-group parallel strategy.
    pub fn resolve_overlaps(&mut self, max_iterations: usize) {
        let self_mutex = std::sync::Arc::new(Mutex::new(self));
        for _ in 0..max_iterations {
            let collision_found = std::sync::Arc::new(AtomicBool::new(false));

            for group_x in 0..3 {
                for group_y in 0..3 {
                    // Collect grid cell keys and their particles for the current group
                    let group_work: Vec<(GridKey, Vec<usize>)> = self_mutex
                        .lock()
                        .unwrap()
                        .grid
                        .grid
                        .iter()
                        .filter(|entry| {
                            let key = *entry.key();
                            key.0.rem_euclid(3) == group_x && key.1.rem_euclid(3) == group_y
                        })
                        .map(|entry| (*entry.key(), entry.value().clone()))
                        .collect();

                    let self_mutex = std::sync::Arc::clone(&self_mutex);
                    let collision_found = std::sync::Arc::clone(&collision_found);
                    group_work
                        .par_iter()
                        .for_each(move |(cell_key, particles_in_cell)| {
                            let mut phx = self_mutex.lock().unwrap();
                            let possible_neighbors = phx.grid.get_possible_neighbors(*cell_key);
                            for &i in particles_in_cell {
                                for &j in &possible_neighbors {
                                    if i <= j {
                                        continue;
                                    }
                                    if phx.resolve_overlap_between(i, j) {
                                        collision_found.store(true, Ordering::Relaxed);
                                    }
                                }
                            }
                        });
                }
            }
            if !Arc::clone(&collision_found).load(Ordering::Relaxed) {
                if !collision_found.load(Ordering::Relaxed) {
                    break;
                }
            }
        }
    }

    /// Apply gravity to all particles.
    pub fn apply_gravity(&mut self) {
        self.particles
            .velocity_y
            .par_iter_mut()
            .for_each(|v| *v += GRAVITY);
    }

    /// Update the simulation state.
    pub fn update(&mut self) {
        self.update_kinematics();
        self.resolve_wall_collisions();
        self.resolve_collisions();
        self.apply_velocity_damping();
        self.apply_gravity();
        self.resolve_overlaps(3);
    }

    pub fn get_drawable_particles(&self) -> Vec<(f32, f32, f32)> {
        self.particles
            .center_x
            .iter()
            .zip(self.particles.center_y.iter())
            .zip(self.particles.radius.iter())
            .map(|((&x, &y), &r)| (x, y, r))
            .collect()
    }
}
