use std::collections::HashMap;
use std::time::{Duration, Instant};

use crate::particles::{ParticleID, Particles};
use crate::sweep::SweepAndPrune;

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
        self.pcls.integrate_v();
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
        self.pcls.apply_gravity();
        self.pcls.apply_velocity_damping();
        self.pcls.integrate_v();

        self.broadphase.update(&self.pcls);
        self.pcls.resolve_wall_collisions();

        self.resolve_overlaps(3);
        self.resolve_collisions();
    }

    pub fn get_drawable_particles(&self) -> (&[f32], &[f32], &[f32]) {
        (&self.pcls.x, &self.pcls.y, &self.pcls.radius)
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
