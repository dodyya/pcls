use std::thread;

use rayon::prelude::*;

use crate::grid::HashGrid;
use crate::particles::ParticleID;

use crate::particles::Particles;
pub struct Phx {
    pub pcls: Particles,
    pub grid: HashGrid,
}
const DT: f32 = 1.0 / 12.0;

impl Phx {
    pub fn new_2k(cell_size: f32) -> Self {
        let mut pcls = Particles::new(30_000);
        pcls.add_10k(0.0, 0.1, cell_size / 2.0, 1.0);
        pcls.add_10k(0.0, 0.0, cell_size / 2.0, 1.0);
        pcls.add_10k(0.0, -0.1, cell_size / 2.0, 1.0);
        let mut grid = HashGrid::new(cell_size);
        for i in 0..pcls.count {
            grid.try_insert(i, pcls.x[i], pcls.y[i]);
        }
        Self { pcls, grid }
    }

    pub fn new(cell_size: f32) -> Self {
        Self {
            pcls: Particles::new(1),
            grid: HashGrid::new(cell_size),
        }
    }

    #[inline(never)]
    pub fn resolve_overlaps(&mut self) {
        thread::spawn(move || {
            self.resolve_some_overlaps(0..self.grid.cell_count);
        });
    }

    pub fn resolve_some_overlaps(&mut self, r: std::ops::Range<usize>) {
        let c = self.grid.cell_count;
        let mut outer: Vec<&[Option<usize>]> = Vec::with_capacity(5);
        for i in r {
            for j in 0..c {
                let inner_ids = &self.grid.map[(i, j)];
                outer.clear();
                outer.push(&self.grid.map[(i, j)]);
                if i < c - 1 {
                    outer.push(&self.grid.map[(i + 1, j)]);
                }
                if j < c - 1 {
                    outer.push(&self.grid.map[(i, j + 1)]);
                }
                if i < c - 1 && j < c - 1 {
                    outer.push(&self.grid.map[(i + 1, j + 1)]);
                }
                if i > 0 && j < c - 1 {
                    outer.push(&self.grid.map[(i - 1, j + 1)]);
                }

                for &in_id in inner_ids.iter().take_while(|x| x.is_some()) {
                    for &out_ids in &outer {
                        for &out_id in out_ids.iter().take_while(|x| x.is_some()) {
                            if let (Some(in_val), Some(out_val)) = (in_id, out_id) {
                                if in_val != out_val {
                                    unsafe {
                                        self.pcls.overlap(in_val, out_val);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[inline(never)]
    pub fn step(&mut self) {
        for _ in 0..10 {
            self.pcls.apply_gravity();
            self.pcls.constrain();
            self.grid.update(&self.pcls);
            self.resolve_overlaps();
            self.pcls.verlet(DT / 10.0);
        }
    }

    pub fn get_drawable_particles(&self) -> (&[f32], &[f32], &[f32]) {
        (&self.pcls.x, &self.pcls.y, &self.pcls.radius)
    }

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        self.pcls.push((x, y, radius, vx, vy, mass));
        self.grid.try_insert(self.pcls.count - 1, x, y);
    }

    pub fn clear(&mut self) {
        self.pcls.clear();
        self.grid.map.clear();
    }
    pub fn toggle_gravity(&mut self) {
        self.pcls.g_toward_center = !self.pcls.g_toward_center;
    }

    pub fn stop(&mut self) {
        self.pcls.stop();
    }
}
