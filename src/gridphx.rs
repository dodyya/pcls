use std::thread;

use rayon::prelude::*;

use crate::grid::Grid;
use crate::particles::ParticleID;

use crate::particles::Particles;
#[derive(Debug)]
pub struct Phx {
    pub pcls: Particles,
    pub grid: Grid,
}
unsafe impl Send for Phx {}
unsafe impl Sync for Phx {}
const DT: f32 = 1.0 / 12.0;

impl Phx {
    pub fn new_2k(cell_size: f32) -> Self {
        let mut pcls = Particles::new(30_000);
        pcls.add_10k(0.0, 0.1, cell_size / 2.0, 1.0);
        pcls.add_10k(0.0, 0.0, cell_size / 2.0, 1.0);
        pcls.add_10k(0.0, -0.1, cell_size / 2.0, 1.0);
        let mut grid = Grid::new(cell_size);
        for i in 0..pcls.count {
            grid.try_insert(i, pcls.x[i], pcls.y[i]);
        }
        Self { pcls, grid }
    }

    pub fn new(cell_size: f32) -> Self {
        Self {
            pcls: Particles::new(1),
            grid: Grid::new(cell_size),
        }
    }

    #[inline(never)]
    pub fn resolve_overlaps(&mut self) {
        let c = self.grid.cell_count;
        // Self::resolve_some_overlaps(c);
    }

    // pub fn resolve_some_overlaps(c: usize, grid: &[usize]) {
    //     let mut outer: Vec<&[Option<usize>]> = Vec::with_capacity(5);
    //     for i in 0..c {
    //         for j in 0..c {
    //             let inner_ids = cell(grid,i, j,];
    //             outer.clear();
    //             outer.push(&grid[ind(i, j)]);
    //             if i < c - 1 {
    //                 outer.push(&grid[ind(i + 1, j)]);
    //             }
    //             if j < c - 1 {
    //                 outer.push(&grid[ind(i, j + 1)]);
    //             }
    //             if i < c - 1 && j < c - 1 {
    //                 outer.push(&grid[ind(i + 1, j + 1)]);
    //             }
    //             if i > 0 && j < c - 1 {
    //                 outer.push(&grid[ind(i - 1, j + 1)]);
    //             }

    //             for &in_id in inner_ids.iter().take_while(|x| x.is_some()) {
    //                 for &out_ids in &outer {
    //                     for &out_id in out_ids.iter().take_while(|x| x.is_some()) {
    //                         if let (Some(in_val), Some(out_val)) = (in_id, out_id) {
    //                             if in_val != out_val {
    //                                 unsafe {
    //                                     self.pcls.overlap(in_val, out_val);
    //                                 }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

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

pub fn cell(subgrid: &[usize], i: usize, j: usize, width: usize, depth: usize) -> &[usize] {
    let start = (j * width + i) * depth;
    &subgrid[start..start + depth]
}
