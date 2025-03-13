use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct Particles {
    pub center_x: Vec<f32>,
    pub center_y: Vec<f32>,
    pub radius: Vec<f32>,
    pub velocity_x: Vec<f32>,
    pub velocity_y: Vec<f32>,
    pub mass: Vec<f32>,
    pub count: usize,
}

impl Particles {
    pub fn new(capacity: usize) -> Self {
        Self {
            center_x: Vec::with_capacity(capacity),
            center_y: Vec::with_capacity(capacity),
            radius: Vec::with_capacity(capacity),
            velocity_x: Vec::with_capacity(capacity),
            velocity_y: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            count: 0,
        }
    }

    pub fn from_particles(particles: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let count = particles.len();
        let mut result = Self::new(count);

        for (cx, cy, r, vx, vy, m) in particles {
            result.center_x.push(cx);
            result.center_y.push(cy);
            result.radius.push(r);
            result.velocity_x.push(vx);
            result.velocity_y.push(vy);
            result.mass.push(m);
        }
        result.count = count;
        result
    }

    fn overlap(&self, i: usize, j: usize) -> (bool, f32) {
        let delta_x = self.center_x[i] - self.center_x[j];
        let delta_y = self.center_y[i] - self.center_y[j];
        let distance_sq = delta_x * delta_x + delta_y * delta_y;
        let radius_sum = self.radius[i] + self.radius[j];
        (distance_sq < radius_sum * radius_sum, distance_sq)
    }

    pub fn clear(&mut self) {
        self.center_x.clear();
        self.center_y.clear();
        self.radius.clear();
        self.velocity_x.clear();
        self.velocity_y.clear();
        self.mass.clear();
        self.count = 0;
    }
}

// Constants
const WALL_DAMPING: f32 = 0.9;
const VELOCITY_DAMPING: f32 = 1.00;
const RESTITUTION: f32 = 0.5;
const GRAVITY_X: f32 = 0.01;
const GRAVITY_Y: f32 = 0.01;
const COLLISION_DAMPING: f32 = 1.00;

pub type GridKey = (i32, i32);

/// The grid stores indices (into the Particles arrays)
pub struct HashGrid {
    cell_size: f32,
    cell_count: i32,
    grid: HashMap<GridKey, Vec<usize>>,
}

impl HashGrid {
    pub fn new(cell_size: f32) -> Self {
        Self {
            cell_size,
            cell_count: (2.0 / cell_size).ceil() as i32,
            grid: HashMap::new(),
        }
    }

    pub fn get_grid_key(&self, x: f32, y: f32) -> GridKey {
        (
            (x / 2.0 * self.cell_count as f32).floor() as i32,
            (y / 2.0 * self.cell_count as f32).floor() as i32,
        )
    }

    /// Inserts a particle index into the cell determined by position.
    pub fn insert(&mut self, index: usize, x: f32, y: f32) {
        let grid_key = self.get_grid_key(x, y);
        self.grid.entry(grid_key).or_default().push(index);
    }

    /// Removes a given particle index from the grid cell corresponding to a given particle position.
    pub fn remove(&mut self, index: usize, x: f32, y: f32) {
        let grid_key = self.get_grid_key(x, y);
        if let Some(vec) = self.grid.get_mut(&grid_key) {
            if let Some(pos) = vec.iter().position(|&i| i == index) {
                vec.swap_remove(pos);
                if vec.is_empty() {
                   self.grid.remove(&grid_key);
                }
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

    /// Update grid membership for a single particle. If its grid key has changed,
    /// remove it from its old cell and reinsert it into the new cell.
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
            // Remove from old cell and insert into new one.
            if let Some(cell) = self.grid.get_mut(&old_key) {
                if let Some(pos) = cell.iter().position(|&i| i == index) {
                    cell.swap_remove(pos);
                    if cell.is_empty() {
                        self.grid.remove(&old_key);
                    }
                }
            }
            self.grid.entry(new_key).or_default().push(index);
        }
    }

    pub fn clear(&mut self) {
        self.grid.clear();
    }
}

/// Central simulation structure: a particles store and a hash grid.
pub struct Phx {
    pub particles: Particles,
    pub grid: HashGrid,
}

impl Phx {
    /// Create a new simulation with a given grid cell size.
    pub fn new(cell_size: f32, particles_data: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let particles = Particles::from_particles(particles_data);
        let mut grid = HashGrid::new(cell_size);
        // Insert all particle indices into the grid.
        for i in 0..particles.count {
            grid.insert(i, particles.center_x[i], particles.center_y[i]);
        }
        Self { particles, grid }
    }

    /// Update kinematics: update positions and grid membership.
    pub fn update_kinematics(&mut self) {
        for i in 0..self.particles.count {
            let old_x = self.particles.center_x[i];
            let old_y = self.particles.center_y[i];

            // Update position based on velocity
            self.particles.center_x[i] += self.particles.velocity_x[i];
            self.particles.center_y[i] += self.particles.velocity_y[i];

            // Update the grid cell if needed.
            self.grid.update_particle(
                i,
                old_x,
                old_y,
                self.particles.center_x[i],
                self.particles.center_y[i],
            );
        }
    }

    /// Resolve collisions with walls.
    pub fn resolve_wall_collisions(&mut self) {
        for i in 0..self.particles.count {
            let old_x = self.particles.center_x[i];
            let old_y = self.particles.center_y[i];
            let radius = self.particles.radius[i];

            // Left wall
            if self.particles.center_x[i] - radius < -1.0 {
                self.particles.center_x[i] = -1.0 + radius;
                self.particles.velocity_x[i] = WALL_DAMPING * self.particles.velocity_x[i].abs();
            }
            // Right wall
            if self.particles.center_x[i] + radius > 1.0 {
                self.particles.center_x[i] = 1.0 - radius;
                self.particles.velocity_x[i] = -WALL_DAMPING * self.particles.velocity_x[i].abs();
            }
            // Bottom wall
            if self.particles.center_y[i] - radius < -1.0 {
                self.particles.center_y[i] = -1.0 + radius;
                self.particles.velocity_y[i] = WALL_DAMPING * self.particles.velocity_y[i].abs();
            }
            // Top wall
            if self.particles.center_y[i] + radius > 1.0 {
                self.particles.center_y[i] = 1.0 - radius;
                self.particles.velocity_y[i] = -WALL_DAMPING * self.particles.velocity_y[i].abs();
            }

            // Update grid membership if needed.
            self.grid.update_particle(
                i,
                old_x,
                old_y,
                self.particles.center_x[i],
                self.particles.center_y[i],
            );
        }
    }

    /// Resolve collisions between particles.
    pub fn resolve_collisions(&mut self) {
        // Get grid keys to iterate over.
        let keys: Vec<GridKey> = self.grid.grid.keys().cloned().collect();
        for cell_key in keys {
            let neighbor_indices = self.grid.get_possible_neighbors(cell_key);
            // For each particle in the current cell, check collisions with all neighbors.
            if let Some(cell_indices) = self.grid.grid.get(&cell_key).cloned() {
                for &i in &cell_indices {
                    // For each neighbor index (avoid self‑collision)
                    for &j in &neighbor_indices {
                        if i <= j {
                            continue;
                        }
                        // Resolve collision between particles[i] and particles[j]
                        let (overlap, _) = self.particles.overlap(i, j);
                        if overlap {
                            self.resolve_collision_between(i, j);
                        }
                    }
                }
            }
        }
    }

    fn resolve_collision_between(&mut self, i: usize, j: usize) {
        // Calculate distance between particles
        let delta_x = self.particles.center_x[i] - self.particles.center_x[j];
        let delta_y = self.particles.center_y[i] - self.particles.center_y[j];
        let distance_sq = delta_x * delta_x + delta_y * delta_y;
        let distance = distance_sq.sqrt();

        // Calculate normal vector
        let (normal_x, normal_y) = if distance != 0.0 {
            (delta_x / distance, delta_y / distance)
        } else {
            (1.0, 0.0)
        };

        // Relative velocity along the normal.
        let relative_velocity_x = self.particles.velocity_x[i] - self.particles.velocity_x[j];
        let relative_velocity_y = self.particles.velocity_y[i] - self.particles.velocity_y[j];
        let velocity_along_normal = relative_velocity_x * normal_x + relative_velocity_y * normal_y;

        if velocity_along_normal > 0.0 {
            return;
        }

        let a_inv = 1.0 / self.particles.mass[i];
        let b_inv = 1.0 / self.particles.mass[j];
        let impulse_mag = -(1.0 + RESTITUTION) * velocity_along_normal / (a_inv + b_inv);
        let impulse_x = normal_x * impulse_mag * COLLISION_DAMPING;
        let impulse_y = normal_y * impulse_mag * COLLISION_DAMPING;

        self.particles.velocity_x[i] += impulse_x * a_inv;
        self.particles.velocity_y[i] += impulse_y * a_inv;
        self.particles.velocity_x[j] -= impulse_x * b_inv;
        self.particles.velocity_y[j] -= impulse_y * b_inv;
    }

    /// Apply velocity damping to all particles.
    pub fn apply_velocity_damping(&mut self) {
        for i in 0..self.particles.count {
            self.particles.velocity_x[i] *= VELOCITY_DAMPING;
            self.particles.velocity_y[i] *= VELOCITY_DAMPING;
        }
    }

    fn resolve_overlap_between(&mut self, i: usize, j: usize) -> bool {
        // Avoid self-collision.
        if i == j {
            return false;
        }

        let (overlap, distance_sq) = self.particles.overlap(i, j);
        if !overlap {
            return false;
        }

        let distance = distance_sq.sqrt();
        let delta_x = self.particles.center_x[i] - self.particles.center_x[j];
        let delta_y = self.particles.center_y[i] - self.particles.center_y[j];

        let (normal_x, normal_y) = if distance != 0.0 {
            (delta_x / distance, delta_y / distance)
        } else {
            (1.0, 0.0)
        };

        let overlap_distance = self.particles.radius[i] + self.particles.radius[j] - distance;
        let correction_x = normal_x * (overlap_distance * 0.5);
        let correction_y = normal_y * (overlap_distance * 0.5);

        // Store old positions for grid membership update.
        let old_pos_i_x = self.particles.center_x[i];
        let old_pos_i_y = self.particles.center_y[i];
        let old_pos_j_x = self.particles.center_x[j];
        let old_pos_j_y = self.particles.center_y[j];

        // Adjust positions.
        self.particles.center_x[i] += correction_x;
        self.particles.center_y[i] += correction_y;
        self.particles.center_x[j] -= correction_x;
        self.particles.center_y[j] -= correction_y;

        // Update grid membership if needed.
        self.grid.update_particle(
            i,
            old_pos_i_x,
            old_pos_i_y,
            self.particles.center_x[i],
            self.particles.center_y[i],
        );
        self.grid.update_particle(
            j,
            old_pos_j_x,
            old_pos_j_y,
            self.particles.center_x[j],
            self.particles.center_y[j],
        );

        true
    }

    /// Iteratively resolve overlaps until no collisions are detected or until max_iterations is reached.
    pub fn resolve_overlaps(&mut self, max_iterations: usize) {
        for _ in 0..max_iterations {
            let mut collision_found = false;
            // Clone grid keys to avoid borrow issues.
            let cell_keys: Vec<GridKey> = self.grid.grid.keys().cloned().collect();
            for key in cell_keys {
                // Clone indices from this cell.
                let cell_indices = self.grid.grid.get(&key).cloned().unwrap_or_default();
                for &i in &cell_indices {
                    let key_i = self
                        .grid
                        .get_grid_key(self.particles.center_x[i], self.particles.center_y[i]);
                    // Get all potential neighbors.
                    let neighbors = self.grid.get_possible_neighbors(key_i);
                    // For each neighbor index, try to resolve a collision.
                    for j in neighbors {
                        if self.resolve_overlap_between(i, j) {
                            collision_found = true;
                        }
                    }
                }
            }
            if !collision_found {
                break;
            }
        }
    }

    
    /// Apply gravity to all particles.
    pub fn apply_gravity(&mut self) {
        for i in 0..self.particles.count {
            //Gravity towards center of screen depending on whether position is positive
            self.particles.velocity_x[i] -= GRAVITY_X * self.particles.center_x[i];
            self.particles.velocity_y[i] -= GRAVITY_Y * self.particles.center_y[i];
        }

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
        let mut result = Vec::with_capacity(self.particles.count);
        for i in 0..self.particles.count {
            result.push((
                self.particles.center_x[i],
                self.particles.center_y[i],
                self.particles.radius[i],
            ));
        }
        result
    }
}
