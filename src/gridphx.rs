use std::collections::HashMap;

use crate::particles::Particles;

// Constants
const WALL_DAMPING: f32 = 0.9;
const VELOCITY_DAMPING: f32 = 1.00;
const RESTITUTION: f32 = 0.5;
const GRAVITY: f32 = 0.01;
const COLLISION_DAMPING: f32 = 1.00;
const GRAVITY_TOWARDS_CENTER: bool = false; // Set to false for downward gravity

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

    /// Returns a vector of indices from neighboring cells that haven't been checked yet.
    /// Uses a spatial optimization to only check cells in the positive direction
    /// (right, down, and down-right diagonal), avoiding duplicate checks.
    pub fn get_possible_neighbors(&self, key: GridKey) -> Vec<usize> {
        let mut neighbors = vec![];

        // Get particles from the current cell
        if let Some(indices) = self.grid.get(&key) {
            neighbors.extend(indices.iter().cloned());
        }

        // Only check these relative positions: right, down-right diagonal, down
        // This ensures we only check each cell pair once
        let relative_positions = [(1, 0), (1, 1), (0, 1), (-1, 1)];

        for (dx, dy) in relative_positions.iter() {
            let cell_key = (key.0 + dx, key.1 + dy);
            if let Some(indices) = self.grid.get(&cell_key) {
                neighbors.extend(indices.iter().cloned());
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
            // Get particles from this cell and neighboring cells in positive directions only
            let all_relevant_indices = self.grid.get_possible_neighbors(cell_key);

            // For each particle in the current cell
            if let Some(current_cell_indices) = self.grid.grid.get(&cell_key).cloned() {
                for (idx, &i) in current_cell_indices.iter().enumerate() {
                    // Check against other particles from same cell (but only in one direction)
                    // This avoids checking (i,j) and (j,i)
                    for &j in current_cell_indices.iter().skip(idx + 1) {
                        // Resolve collision between particles[i] and particles[j]
                        let (overlap, _) = self.particles.overlap(i, j);
                        if overlap {
                            self.resolve_collision_between(i, j);
                        }
                    }

                    // Check against particles from neighboring cells
                    // We only need to check against particles from other cells
                    for &j in &all_relevant_indices {
                        // Skip particles from the same cell - we already checked them above
                        if current_cell_indices.contains(&j) {
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
                // Get particles from this cell and neighboring cells in positive directions only
                let all_relevant_indices = self.grid.get_possible_neighbors(key);

                // For each particle in the current cell
                if let Some(current_cell_indices) = self.grid.grid.get(&key).cloned() {
                    for (idx, &i) in current_cell_indices.iter().enumerate() {
                        // Check against other particles from same cell (but only in one direction)
                        // This avoids checking (i,j) and (j,i)
                        for &j in current_cell_indices.iter().skip(idx + 1) {
                            if self.resolve_overlap_between(i, j) {
                                collision_found = true;
                            }
                        }

                        // Check against particles from neighboring cells
                        // We only need to check against particles from other cells
                        for &j in &all_relevant_indices {
                            // Skip particles from the same cell - we already checked them above
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

    /// Apply gravity to all particles.
    pub fn apply_gravity(&mut self) {
        for i in 0..self.particles.count {
            if GRAVITY_TOWARDS_CENTER {
                // Gravity towards center of screen
                self.particles.velocity_x[i] -= GRAVITY * self.particles.center_x[i];
                self.particles.velocity_y[i] -= GRAVITY * self.particles.center_y[i];
            } else {
                // Downward gravity
                self.particles.velocity_y[i] -= 0.2 * GRAVITY;
            }
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

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        self.particles.push((x, y, radius, vx, vy, mass));
        self.grid.insert(self.particles.count - 1, x, y);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use std::time::{Duration, Instant};

    fn create_test_simulation(particle_count: usize) -> Phx {
        let mut rng = rand::thread_rng();
        let mut particles = Vec::with_capacity(particle_count);

        for _ in 0..particle_count {
            let x = rng.gen_range(-0.9..0.9);
            let y = rng.gen_range(-0.9..0.9);
            let radius = rng.gen_range(0.01..0.05);
            let vx = rng.gen_range(-0.01..0.01);
            let vy = rng.gen_range(-0.01..0.01);
            let mass = radius * radius * std::f32::consts::PI;

            particles.push((x, y, radius, vx, vy, mass));
        }

        Phx::new(0.1, particles)
    }

    fn time_function<F>(name: &str, mut f: F) -> Duration
    where
        F: FnMut(),
    {
        let start = Instant::now();
        f();
        let duration = start.elapsed();
        println!("{} took: {:?}", name, duration);
        duration
    }

    #[test]
    fn benchmark_simulation_steps() {
        let particle_counts = [100, 500, 1000, 5000];
        let mut results = HashMap::new();

        for &count in &particle_counts {
            println!("\nBenchmarking with {} particles", count);
            let mut sim = create_test_simulation(count);

            let mut step_times = HashMap::new();

            // Time each step individually
            step_times.insert(
                "update_kinematics",
                time_function("update_kinematics", || sim.update_kinematics()),
            );
            step_times.insert(
                "resolve_wall_collisions",
                time_function("resolve_wall_collisions", || sim.resolve_wall_collisions()),
            );
            step_times.insert(
                "resolve_collisions",
                time_function("resolve_collisions", || sim.resolve_collisions()),
            );
            step_times.insert(
                "apply_velocity_damping",
                time_function("apply_velocity_damping", || sim.apply_velocity_damping()),
            );
            step_times.insert(
                "apply_gravity",
                time_function("apply_gravity", || sim.apply_gravity()),
            );
            step_times.insert(
                "resolve_overlaps",
                time_function("resolve_overlaps", || sim.resolve_overlaps(3)),
            );

            // Time full update
            let full_update_time = time_function("full update", || sim.update());
            step_times.insert("full_update", full_update_time);

            results.insert(count, step_times);
        }

        // Print summary
        println!("\n=== PERFORMANCE SUMMARY ===");
        for &count in &particle_counts {
            println!("\nParticle count: {}", count);
            if let Some(times) = results.get(&count) {
                for (step, duration) in times {
                    println!("  {} - {:?}", step, duration);
                }
            }
        }
    }

    #[test]
    fn test_particle_count_until_latency_threshold() {
        let mut rng = rand::thread_rng();
        let mut sim = Phx::new(0.1, vec![]);
        let target_latency = Duration::from_millis(5);
        let mut particle_count = 0;
        let mut duration = Duration::from_millis(0);
        let batch_size = 50; // Add particles in batches

        println!("Adding particles until latency reaches 5ms...");

        while duration < target_latency {
            // Add a batch of particles
            for _ in 0..batch_size {
                let x = 0.0;
                let y = 0.0;
                let radius = rng.gen_range(0.01..0.05);
                let vx = rng.gen_range(-0.01..0.01);
                let vy = rng.gen_range(-0.01..0.01);
                let mass = radius * radius * std::f32::consts::PI;
                sim.add_particle(x, y, radius, vx, vy, mass);
                particle_count += 1;
            }

            // Measure the time it takes to update the simulation
            let start = Instant::now();
            sim.update();
            duration = start.elapsed();

            // Print progress every 500 particles
            if particle_count % 500 == 0 {
                println!(
                    "Particles: {}, Current latency: {:?}",
                    particle_count, duration
                );
            }
        }

        println!("\n=== LATENCY THRESHOLD REACHED ===");
        println!("Number of particles: {}", particle_count);
        println!("Final update latency: {:?}", duration);

        // Assert that we've reached the threshold
        assert!(duration >= target_latency);
    }
}
