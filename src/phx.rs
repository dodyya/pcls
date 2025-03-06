use std::{
    collections::HashMap,
    ops::{Add, AddAssign, Div, Mul, Sub, SubAssign},
};

#[derive(Debug, Clone, PartialEq)]
pub struct Particle {
    pub center: Vec2<f32>,
    pub radius: f32,
    pub velocity: Vec2<f32>,
    pub mass: f32,
}

impl Particle {
    fn overlap(&self, other: &Particle) -> (bool, f32) {
        let delta = self.center - other.center;
        let distance_sq = delta.magnitude_sq();
        let radius_sum = self.radius + other.radius;
        (distance_sq < radius_sum * radius_sum, distance_sq)
    }
}

// Constant wall damping coefficient:
const WALL_DAMPING: f32 = 0.9;
const VELOCITY_DAMPING: f32 = 0.99;
const RESTITUTION: f32 = 0.5;
const GRAVITY: Vec2<f32> = Vec2 { x: 0.0, y: -0.001 };
const COLLISION_DAMPING: f32 = 1.00;

pub type GridKey = Vec2<i32>;

/// The grid now stores indices (into the Phx.particles Vec)
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
        Vec2 {
            x: (x / 2.0 * self.cell_count as f32).floor() as i32,
            y: (y / 2.0 * self.cell_count as f32).floor() as i32,
        }
    }

    /// Inserts a particle index into the cell determined by p’s position.
    pub fn insert(&mut self, index: usize, p: &Particle) {
        let grid_key = self.get_grid_key(p.center.x, p.center.y);
        self.grid.entry(grid_key).or_default().push(index);
    }

    /// Removes a given particle index from the grid cell corresponding to a given particle position.
    pub fn remove(&mut self, index: usize, p: &Particle) {
        let grid_key = self.get_grid_key(p.center.x, p.center.y);
        if let Some(vec) = self.grid.get_mut(&grid_key) {
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
                let cell_key = key + Vec2 { x: i, y: j };
                if let Some(indices) = self.grid.get(&cell_key) {
                    neighbors.extend(indices.iter().cloned());
                }
            }
        }
        neighbors
    }

    /// Update grid membership for a single particle. If its grid key has changed,
    /// remove it from its old cell and reinsert it into the new cell.
    pub fn update_particle(&mut self, index: usize, old_pos: &Vec2<f32>, new_particle: &Particle) {
        let old_key = self.get_grid_key(old_pos.x, old_pos.y);
        let new_key = self.get_grid_key(new_particle.center.x, new_particle.center.y);
        if old_key != new_key {
            // Remove from old cell and insert into new one.
            if let Some(cell) = self.grid.get_mut(&old_key) {
                if let Some(pos) = cell.iter().position(|&i| i == index) {
                    cell.swap_remove(pos);
                }
            }
            self.grid.entry(new_key).or_default().push(index);
        }
    }
}

/// Central simulation structure: a particle store and a hash grid.
pub struct Phx {
    pub particles: Vec<Particle>,
    pub grid: HashGrid,
}

impl Phx {
    /// Create a new simulation with a given grid cell size.
    pub fn new(cell_size: f32, particles: Vec<Particle>) -> Self {
        let mut grid = HashGrid::new(cell_size);
        // Insert all particle indices into the grid.
        for (i, p) in particles.iter().enumerate() {
            grid.insert(i, p);
        }
        Self { particles, grid }
    }

    /// Update kinematics: update positions and grid membership.
    pub fn update_kinematics(&mut self) {
        for (i, p) in self.particles.iter_mut().enumerate() {
            let old_pos = p.center;
            p.center += p.velocity;
            // Update the grid cell if needed.
            self.grid.update_particle(i, &old_pos, p);
        }
    }

    /// Resolve collisions with walls.
    pub fn resolve_wall_collisions(&mut self) {
        // We determine affected particles by checking the particle position.
        for (i, p) in self.particles.iter_mut().enumerate() {
            let old_pos = p.center;
            if p.center.x - p.radius < -1.0 {
                p.center.x = -1.0 + p.radius;
                p.velocity.x = WALL_DAMPING * p.velocity.x.abs();
            }
            if p.center.x + p.radius > 1.0 {
                p.center.x = 1.0 - p.radius;
                p.velocity.x = -WALL_DAMPING * p.velocity.x.abs();
            }
            if p.center.y - p.radius < -1.0 {
                p.center.y = -1.0 + p.radius;
                p.velocity.y = WALL_DAMPING * p.velocity.y.abs();
            }
            if p.center.y + p.radius > 1.0 {
                p.center.y = 1.0 - p.radius;
                p.velocity.y = -WALL_DAMPING * p.velocity.y.abs();
            }
            // Update grid membership if needed.
            self.grid.update_particle(i, &old_pos, p);
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
                        // (Note: this simple scheme may resolve collisions twice.)
                        let (overlap, _) = self.particles[i].overlap(&self.particles[j]);
                        if overlap {
                            self.resolve_collision_between(i, j);
                        }
                    }
                }
            }
        }
    }

    fn resolve_collision_between(&mut self, i: usize, j: usize) {
        // println!("Resolving collision between particles {} and {}", i, j);
        let distance = (self.particles[i].center - self.particles[j].center)
            .magnitude_sq()
            .sqrt();
        let normal = if distance != 0.0 {
            (self.particles[i].center - self.particles[j].center) * (1.0 / distance)
        } else {
            Vec2 { x: 1.0, y: 0.0 }
        };
        // Relative velocity along the normal.
        let relative_velocity = self.particles[i].velocity - self.particles[j].velocity;
        let velocity_along_normal = relative_velocity.dot(normal);
        if velocity_along_normal > 0.0 {
            return;
        }
        let a_inv = 1.0 / self.particles[i].mass;
        let b_inv = 1.0 / self.particles[j].mass;
        let impulse_mag = -(1.0 + RESTITUTION) * velocity_along_normal / (a_inv + b_inv);
        let impulse = normal * impulse_mag * COLLISION_DAMPING;
        self.particles[i].velocity = self.particles[i].velocity + impulse * a_inv;
        self.particles[j].velocity = self.particles[j].velocity - impulse * b_inv;
    }

    /// Apply velocity damping to all particles.
    pub fn apply_velocity_damping(&mut self) {
        for p in self.particles.iter_mut() {
            p.velocity = p.velocity * VELOCITY_DAMPING;
        }
    }

    fn resolve_overlap_between(&mut self, i: usize, j: usize) -> bool {
        // Avoid self-collision.
        if i == j {
            return false;
        }
        let (overlap, distance_sq) = self.particles[i].overlap(&self.particles[j]);
        if !overlap {
            return false;
        }
        let distance = distance_sq.sqrt();
        let normal = if distance != 0.0 {
            (self.particles[i].center - self.particles[j].center) * (1.0 / distance)
        } else {
            Vec2 { x: 1.0, y: 0.0 }
        };
        let overlap_distance = self.particles[i].radius + self.particles[j].radius - distance;
        let correction = normal * (overlap_distance * 0.5);

        // Store old positions for grid membership update.
        let old_pos_i = self.particles[i].center;
        let old_pos_j = self.particles[j].center;

        // Adjust positions.
        self.particles[i].center += correction;
        self.particles[j].center -= correction;

        // Update grid membership if needed.
        self.grid.update_particle(i, &old_pos_i, &self.particles[i]);
        self.grid.update_particle(j, &old_pos_j, &self.particles[j]);
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
                        .get_grid_key(self.particles[i].center.x, self.particles[i].center.y);
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
        for p in self.particles.iter_mut() {
            p.velocity += GRAVITY;
        }
    }

    /// Update the simulation state.
    pub fn update(&mut self) {
        self.update_kinematics();
        self.resolve_wall_collisions();
        self.resolve_collisions();
        self.apply_velocity_damping();
        self.apply_gravity();
        self.resolve_overlaps(20);
    }

    pub fn get_drawable_particles(&self) -> Vec<(f32, f32, f32)> {
        return self
            .particles
            .iter()
            .map(|p| (p.center.x, p.center.y, p.radius))
            .collect();
    }
}

// ===== Vec2 Implementation =====

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Vec2<T: Add<Output = T> + Copy> {
    pub x: T,
    pub y: T,
}

impl<T: Add<Output = T> + Copy> Add for Vec2<T> {
    type Output = Vec2<T>;

    fn add(self, other: Vec2<T>) -> Vec2<T> {
        Vec2 {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl<T: Sub<Output = T> + Add<Output = T> + Copy> Sub for Vec2<T> {
    type Output = Vec2<T>;

    fn sub(self, other: Vec2<T>) -> Vec2<T> {
        Vec2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl<T: Add<Output = T> + Sub<Output = T> + Copy> SubAssign for Vec2<T> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}
impl<T: Add<Output = T> + Copy> AddAssign for Vec2<T> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<T: Mul<Output = T> + Copy + Add<Output = T>> Vec2<T> {
    fn magnitude_sq(&self) -> T {
        self.x * self.x + self.y * self.y
    }

    fn dot(&self, other: Vec2<T>) -> T {
        self.x * other.x + self.y * other.y
    }
}

impl<T> Mul<T> for Vec2<T>
where
    T: Add<Output = T> + Mul<Output = T> + Copy,
{
    type Output = Self;

    fn mul(self, scalar: T) -> Self::Output {
        Self {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hashgrid_new() {
        let grid = HashGrid::new(10.0);
        assert_eq!(grid.cell_size, 10.0);
        assert!(grid.grid.is_empty());
    }

    #[test]
    fn test_hashgrid_get_grid_key() {
        let grid = HashGrid::new(0.1);
        let key = grid.get_grid_key(0.1, 0.2);
        assert_eq!(key, Vec2 { x: 1, y: 2 });
    }

    #[test]
    fn test_phx_integration() {
        // Create two particles.
        let particle1 = Particle {
            center: Vec2 { x: 15.0, y: 25.0 },
            radius: 1.0,
            velocity: Vec2 { x: 1.0, y: 1.0 },
            mass: 1.0,
        };
        let particle2 = Particle {
            center: Vec2 { x: 25.0, y: 35.0 },
            radius: 1.0,
            velocity: Vec2 { x: 1.0, y: 1.0 },
            mass: 1.0,
        };

        let particles = vec![particle1.clone(), particle2.clone()];
        let mut phx = Phx::new(10.0, particles);

        // Check that grid contains both particle indices.
        let key = phx.grid.get_grid_key(15.0, 25.0);
        let neighbors = phx.grid.get_possible_neighbors(key);
        assert!(neighbors.contains(&0));
        // (Depending on positions, particle2 might be in a different cell.)
    }
}
