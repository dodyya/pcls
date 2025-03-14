use std::collections::HashMap;
use std::time::{Duration, Instant};

use crate::particles::Particles;

// Constants
const WALL_DAMPING: f32 = 0.9;
const VELOCITY_DAMPING: f32 = 1.00;
const RESTITUTION: f32 = 0.5;
const GRAVITY: f32 = 0.01;
const COLLISION_DAMPING: f32 = 1.00;
const GRAVITY_TOWARDS_CENTER: bool = false; // Set to false for downward gravity

/// Particle index with relevant coordinate for sweep and prune sorting
#[derive(Debug, Clone, Copy)]
struct AxisEntry {
    index: usize,
    min: f32,
    max: f32,
}

/// A sweep-and-prune implementation for efficient broadphase collision detection
pub struct SweepAndPrune {
    /// Entries sorted by x-axis
    entries_x: Vec<AxisEntry>,
    /// Entries sorted by y-axis
    entries_y: Vec<AxisEntry>,
}

impl SweepAndPrune {
    pub fn new(particles: &Particles) -> Self {
        let mut entries_x = Vec::with_capacity(particles.count);
        let mut entries_y = Vec::with_capacity(particles.count);

        for i in 0..particles.count {
            let radius = particles.radius[i];

            // Create entry for x-axis
            entries_x.push(AxisEntry {
                index: i,
                min: particles.center_x[i] - radius,
                max: particles.center_x[i] + radius,
            });

            // Create entry for y-axis
            entries_y.push(AxisEntry {
                index: i,
                min: particles.center_y[i] - radius,
                max: particles.center_y[i] + radius,
            });
        }

        // Initial sort of entries
        entries_x.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        entries_y.sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());

        Self {
            entries_x,
            entries_y,
        }
    }

    /// Clear all entries from the data structure
    pub fn clear(&mut self) {
        self.entries_x.clear();
        self.entries_y.clear();
    }

    /// Update all entries with new particle positions and radii
    pub fn update(&mut self, particles: &Particles) {
        // Update x-axis entries
        for entry in &mut self.entries_x {
            let i = entry.index;
            entry.min = particles.center_x[i] - particles.radius[i];
            entry.max = particles.center_x[i] + particles.radius[i];
        }

        // Update y-axis entries
        for entry in &mut self.entries_y {
            let i = entry.index;
            entry.min = particles.center_y[i] - particles.radius[i];
            entry.max = particles.center_y[i] + particles.radius[i];
        }

        // Re-sort both arrays
        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        self.entries_y
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }

    /// Get potentially colliding pairs by sweeping along the x-axis
    pub fn get_potential_collisions(&self) -> Vec<(usize, usize)> {
        let mut potential_pairs = Vec::new();
        let entries = &self.entries_x; // We'll just use x-axis for simplicity

        // For each entry, find all entries that could potentially overlap
        for i in 0..entries.len() {
            let entry_i = &entries[i];

            // Sweep forward to find all entries that overlap with entry_i
            for j in (i + 1)..entries.len() {
                let entry_j = &entries[j];

                // If the min of j is past the max of i, no more collisions are possible
                if entry_j.min > entry_i.max {
                    break;
                }

                // Potential collision detected
                potential_pairs.push((entry_i.index, entry_j.index));
            }
        }

        potential_pairs
    }

    /// Clear and rebuild the whole data structure
    pub fn rebuild(&mut self, particles: &Particles) {
        self.entries_x.clear();
        self.entries_y.clear();

        for i in 0..particles.count {
            let radius = particles.radius[i];

            self.entries_x.push(AxisEntry {
                index: i,
                min: particles.center_x[i] - radius,
                max: particles.center_x[i] + radius,
            });

            self.entries_y.push(AxisEntry {
                index: i,
                min: particles.center_y[i] - radius,
                max: particles.center_y[i] + radius,
            });
        }

        self.entries_x
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
        self.entries_y
            .sort_by(|a, b| a.min.partial_cmp(&b.min).unwrap());
    }
}

/// Central simulation structure: a particles store and a sweep-and-prune broadphase
pub struct Phx {
    pub particles: Particles,
    pub broadphase: SweepAndPrune,
}

impl Phx {
    /// Create a new simulation
    pub fn new(_cell_size: f32, particles_data: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let particles = Particles::from_particles(particles_data);
        let broadphase = SweepAndPrune::new(&particles);

        Self {
            particles,
            broadphase,
        }
    }

    /// Update kinematics: update positions
    pub fn update_kinematics(&mut self) {
        for i in 0..self.particles.count {
            // Update position based on velocity
            self.particles.center_x[i] += self.particles.velocity_x[i];
            self.particles.center_y[i] += self.particles.velocity_y[i];
        }

        // Update the broadphase with new positions
        self.broadphase.update(&self.particles);
    }

    /// Resolve collisions with walls.
    pub fn resolve_wall_collisions(&mut self) {
        for i in 0..self.particles.count {
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
        }

        // Update the broadphase after wall collisions
        self.broadphase.update(&self.particles);
    }

    /// Resolve collisions between particles.
    pub fn resolve_collisions(&mut self) {
        // Get potentially colliding pairs from the broadphase
        let potential_collisions = self.broadphase.get_potential_collisions();

        // For each pair, check actual collision and resolve it
        for (i, j) in potential_collisions {
            // Resolve collision between particles[i] and particles[j]
            let (overlap, _) = self.particles.overlap(i, j);
            if overlap {
                self.resolve_collision_between(i, j);
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

        // Adjust positions.
        self.particles.center_x[i] += correction_x;
        self.particles.center_y[i] += correction_y;
        self.particles.center_x[j] -= correction_x;
        self.particles.center_y[j] -= correction_y;

        true
    }

    /// Iteratively resolve overlaps until no collisions are detected or until max_iterations is reached.
    pub fn resolve_overlaps(&mut self, max_iterations: usize) {
        for _ in 0..max_iterations {
            let mut collision_found = false;

            // Get potentially colliding pairs from the broadphase
            let potential_collisions = self.broadphase.get_potential_collisions();

            // For each pair, try to resolve overlaps
            for (i, j) in potential_collisions {
                if self.resolve_overlap_between(i, j) {
                    collision_found = true;
                }
            }

            // Update broadphase after position changes
            self.broadphase.update(&self.particles);

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
}
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use std::time::Instant;

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
}
