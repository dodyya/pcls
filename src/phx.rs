use std::sync::Arc;
use std::sync::RwLock;
use std::thread;

use atomic_float::AtomicF32;
use rayon::prelude::*;

use crate::grid::Grid;
use crate::maybe_id::MaybeID;
use crate::particles::O;

use crate::array::Array3D;
use crate::particles::Particles;
#[derive(Debug)]
pub struct Phx {
    pub pcls: Arc<Particles>,
    pub grid: Arc<Grid>,
}
unsafe impl Send for Phx {}
unsafe impl Sync for Phx {}
const DT: f32 = 1.0 / 12.0;
const NUM_ITERS: usize = 10;
const NUM_THREADS: usize = 10;

impl Phx {
    pub fn new_2k(cell_size: f32) -> Self {
        let mut pcls = Particles::new(30_000);
        pcls.add_10k(0.0, 0.1, cell_size / 2.0, 1.0);
        pcls.add_10k(0.0, 0.0, cell_size / 2.0, 1.0);
        pcls.add_10k(0.0, -0.1, cell_size / 2.0, 1.0);
        let grid = Grid::new(cell_size);
        for i in 0..pcls.count {
            grid.try_insert(i, pcls.x[i].load(O), pcls.y[i].load(O));
        }
        Self {
            pcls: Arc::new(pcls),
            grid: Arc::new(grid),
        }
    }

    pub fn new(cell_size: f32) -> Self {
        let grid = Grid::new(cell_size);
        Self {
            pcls: Arc::new(Particles::new(1)),
            grid: Arc::new(grid),
        }
    }

    pub fn resolve_overlaps(&mut self, n_threads: usize) {
        let c = self.grid.cell_count;
        if c % n_threads != 0 {
            panic!("Cells to a side is not divisible by number of threads");
        }

        thread::scope(|s| {
            let thread_width = c / n_threads;
            for n in 0..n_threads {
                let arc_grid = Arc::clone(&self.grid);
                let arc_pcls = Arc::clone(&self.pcls);
                s.spawn(move || {
                    Self::col_overlap(
                        n * thread_width..n * thread_width + thread_width / 2,
                        c,
                        arc_grid,
                        arc_pcls,
                    )
                });
            }
        });

        thread::scope(|s| {
            let thread_width = c / n_threads;
            for n in 0..n_threads {
                let arc_grid = Arc::clone(&self.grid);
                let arc_pcls = Arc::clone(&self.pcls);
                s.spawn(move || {
                    Self::col_overlap(
                        n * thread_width + thread_width / 2..(n + 1) * thread_width,
                        c,
                        arc_grid,
                        arc_pcls,
                    )
                });
            }
        });
    }

    #[inline(never)]
    pub fn step(&mut self) {
        for _ in 0..NUM_ITERS {
            self.pcls.apply_gravity();
            self.pcls.constrain();
            self.grid.update(&self.pcls);
            self.resolve_overlaps(NUM_THREADS);
            self.pcls.verlet(DT / 10.0);
        }
    }

    pub fn get_drawable_particles(&self) -> (&[AtomicF32], &[AtomicF32], &[AtomicF32]) {
        (
            self.pcls.x.as_ref(),
            self.pcls.y.as_ref(),
            self.pcls.r.as_ref(),
        )
    }

    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, vx: f32, vy: f32, mass: f32) {
        Arc::get_mut(&mut self.pcls)
            .unwrap()
            .push((x, y, radius, vx, vy, mass));
        self.grid.try_insert(self.pcls.count - 1, x, y);
    }

    pub fn clear(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().clear();
        self.grid.map.clear();
    }
    pub fn toggle_gravity(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().g_toward_center = !self.pcls.g_toward_center;
    }

    pub fn stop(&mut self) {
        self.pcls.stop();
    }

    fn col_overlap(
        x_range: std::ops::Range<usize>,
        c: usize,
        // n: usize,
        arc_grid: Arc<Grid>,
        arc_pcls: Arc<Particles>,
        // thread_width: usize,
    ) {
        let grid_ref = arc_grid.as_ref();
        let mut outer: Vec<&[MaybeID]> = Vec::with_capacity(4);
        for i in x_range {
            for j in (0..c).rev() {
                let inner_ids = &grid_ref.map[(i, j)];
                outer.clear();
                if i < c - 1 {
                    outer.push(&grid_ref.map[(i + 1, j)]);
                }
                if j < c - 1 {
                    outer.push(&grid_ref.map[(i, j + 1)]);
                }
                if i < c - 1 && j < c - 1 {
                    outer.push(&grid_ref.map[(i + 1, j + 1)]);
                }
                if i > 0 && j < c - 1 {
                    outer.push(&grid_ref.map[(i - 1, j + 1)]);
                }

                for (num, in_id) in inner_ids.iter().take_while(|x| x.is_some()).enumerate() {
                    for &out_ids in &outer {
                        for out_id in out_ids.iter().take_while(|x| x.is_some()) {
                            let (in_val, out_val) = (in_id.id().unwrap(), out_id.id().unwrap());

                            if in_val != out_val {
                                arc_pcls.overlap(in_val, out_val);
                            }
                        }
                    }

                    for other_id in inner_ids.iter().skip(num).take_while(|x| x.is_some()) {
                        let (in_val, out_val) = (in_id.id().unwrap(), other_id.id().unwrap());

                        if in_val != out_val {
                            arc_pcls.overlap(in_val, out_val);
                        }
                    }
                }
            }
        }
    }
}
