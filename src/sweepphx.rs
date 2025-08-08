use std::collections::HashMap;
use std::time::{Duration, Instant};

use crate::particles::{ParticleID, Particles};
use crate::sweep::SweepAndPrune;

const RESTITUTION: f32 = 0.5;
const GRAVITY: f32 = 0.01;
const COLLISION_DAMPING: f32 = 1.00;
const GRAVITY_TOWARDS_CENTER: bool = false;

pub struct Phx {
    pub pcls: Particles,
    pub broadphase: SweepAndPrune,
}

impl Phx {
    pub fn new(_cell_size: f32, particles_data: Vec<(f32, f32, f32, f32, f32, f32)>) -> Self {
        let particles = Particles::from_particles(particles_data);
        let broadphase = SweepAndPrune::new(&particles);

        Self {
            pcls: particles,
            broadphase,
        }
    }

    pub fn update_kinematics(&mut self) {
        self.pcls.update_kinematics();
        self.broadphase.update(&self.pcls);
    }

    pub fn resolve_wall_collisions(&mut self) {
        self.pcls.resolve_wall_collisions();
        self.broadphase.update(&self.pcls);
    }

    pub fn resolve_collisions(&mut self) {
        let potential_collisions = self.broadphase.get_potential_collisions();

        for (i, j) in potential_collisions {
            let (overlap, _) = self.pcls.overlap_info(i, j);
            if overlap {
                self.pcls.collision(i, j);
            }
        }
    }

    fn resolve_collision_between(&mut self, i: usize, j: usize) {
        let delta_x = self.pcls.x[i] - self.pcls.x[j];
        let delta_y = self.pcls.y[i] - self.pcls.y[j];
        let distance_sq = delta_x * delta_x + delta_y * delta_y;
        let distance = distance_sq.sqrt();

        let (normal_x, normal_y) = if distance != 0.0 {
            (delta_x / distance, delta_y / distance)
        } else {
            (1.0, 0.0)
        };

        let relative_velocity_x = self.pcls.vx[i] - self.pcls.vx[j];
        let relative_velocity_y = self.pcls.vy[i] - self.pcls.vy[j];
        let velocity_along_normal = relative_velocity_x * normal_x + relative_velocity_y * normal_y;

        if velocity_along_normal > 0.0 {
            return;
        }

        let a_inv = 1.0 / self.pcls.mass[i];
        let b_inv = 1.0 / self.pcls.mass[j];
        let impulse_mag = -(1.0 + RESTITUTION) * velocity_along_normal / (a_inv + b_inv);
        let impulse_x = normal_x * impulse_mag * COLLISION_DAMPING;
        let impulse_y = normal_y * impulse_mag * COLLISION_DAMPING;

        self.pcls.vx[i] += impulse_x * a_inv;
        self.pcls.vy[i] += impulse_y * a_inv;
        self.pcls.vx[j] -= impulse_x * b_inv;
        self.pcls.vy[j] -= impulse_y * b_inv;
    }

    pub fn resolve_overlaps(&mut self, max_iterations: usize) {
        for _ in 0..max_iterations {
            let mut collision_found = false;

            let potential_collisions = self.broadphase.get_potential_collisions();

            for (i, j) in potential_collisions {
                if self.pcls.overlap(i, j) {
                    collision_found = true;
                }
            }

            self.broadphase.update(&self.pcls);

            if !collision_found {
                break;
            }
        }
    }

    pub fn step(&mut self) {
        self.pcls.update_kinematics(); //Changes position, velocity
        self.pcls.resolve_wall_collisions(); //Changes position, velocity
        self.broadphase.update(&self.pcls);
        self.resolve_collisions(); //Changes velocity
        self.pcls.apply_velocity_damping(); //Changes velocity
        self.pcls.apply_gravity(); //Changes velocity
        self.resolve_overlaps(3); //Changes position
    }

    pub fn get_drawable_particles(&self) -> Vec<(f32, f32, f32)> {
        let mut result = Vec::with_capacity(self.pcls.count);
        for i in 0..self.pcls.count {
            result.push((self.pcls.x[i], self.pcls.y[i], self.pcls.radius[i]));
        }
        result
    }

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        self.pcls.x.push(x);
        self.pcls.y.push(y);
        self.pcls.radius.push(radius);
        self.pcls.vx.push(vx);
        self.pcls.vy.push(vy);
        self.pcls.mass.push(mass);
        self.pcls.count += 1;

        self.broadphase.rebuild(&self.pcls);
    }

    pub fn clear(&mut self) {
        self.pcls.clear();
        self.broadphase.clear();
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
}
