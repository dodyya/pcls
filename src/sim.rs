use crate::grid::Grid;
use crate::maybe_id::MaybeID;
use crate::particles::Particles;
use rand::Rng;
use std::sync::Arc;
use std::thread;

const DT: f32 = 1.0 / 60.0;
const SUBSTEPS: usize = 12;
const NUM_THREADS: usize = 8;
const GRAVITY: f32 = 1.0;
const WASHING_MACHINE: bool = false;
const RESTITUTION: f32 = 1.0;
const ANTI_BLACK_HOLE: f32 = 0.5;
const MAX_V: f32 = 0.01;
const GRID_DEPTH: usize = 3;
const VELOCITY_DAMPING: f32 = 0.99999999;

#[derive(Debug)]
pub struct Simulation {
    pub pcls: Arc<Particles>,
    pub grid: Arc<Grid>,
}

impl Simulation {
    pub fn new(cell_size: f32) -> Self {
        let grid = Grid::new(cell_size, GRID_DEPTH);
        Self {
            pcls: Arc::new(Particles::new(1)),
            grid: Arc::new(grid),
        }
    }
    #[inline(never)]
    pub fn step(&mut self) {
        for _ in 0..SUBSTEPS {
            Self::apply_gravity(&self.pcls);
            Self::constrain(&self.pcls);
            self.grid.update(&self.pcls);
            Self::resolve_overlaps(&self.grid, &self.pcls, NUM_THREADS);
            Self::verlet(&self.pcls, DT / SUBSTEPS as f32);
        }
    }

    #[inline(never)]
    pub fn resolve_overlaps(grid: &Arc<Grid>, pcls: &Arc<Particles>, n_threads: usize) {
        let c = grid.cell_count;
        if c % n_threads != 0 {
            panic!("Cells to a side is not divisible by number of threads");
        }

        thread::scope(|s| {
            let thread_width = c / n_threads;
            for n in 0..n_threads {
                let arc_grid = Arc::clone(grid);
                let arc_pcls = Arc::clone(pcls);
                s.spawn(move || {
                    let range = n * thread_width..n * thread_width + thread_width / 2;
                    Self::overlap_chunk(range, c, arc_grid, arc_pcls)
                });
            }
        });

        thread::scope(|s| {
            let thread_width = c / n_threads;
            for n in 0..n_threads {
                let arc_grid = Arc::clone(grid);
                let arc_pcls = Arc::clone(pcls);
                s.spawn(move || {
                    let range = n * thread_width + thread_width / 2..(n + 1) * thread_width;
                    Self::overlap_chunk(range, c, arc_grid, arc_pcls)
                });
            }
        });
    }

    #[inline(never)]
    fn overlap_chunk(
        x_range: std::ops::Range<usize>,
        c: usize,
        arc_grid: Arc<Grid>,
        arc_pcls: Arc<Particles>,
    ) {
        let grid_ref = arc_grid.as_ref();
        let mut outer: Vec<&[MaybeID]> = Vec::with_capacity(4);
        for i in x_range {
            for j in (0..c).rev() {
                let inner_ids = &grid_ref.map[(i, j)];
                outer.clear();
                if i > 0 && j < c - 1 {
                    outer.push(&grid_ref.map[(i - 1, j + 1)]);
                }
                if i < c - 1 {
                    outer.push(&grid_ref.map[(i + 1, j)]);
                }
                if j < c - 1 {
                    outer.push(&grid_ref.map[(i, j + 1)]);
                }
                if i < c - 1 && j < c - 1 {
                    outer.push(&grid_ref.map[(i + 1, j + 1)]);
                }

                for (num, in_id) in inner_ids.iter().enumerate() {
                    if in_id.is_none() {
                        break;
                    }
                    for &out_ids in &outer {
                        for out_id in out_ids {
                            if out_id.is_none() {
                                break;
                            }
                            let (in_val, out_val) = (in_id.id().unwrap(), out_id.id().unwrap());
                            if in_val != out_val {
                                Self::overlap(arc_pcls.as_ref(), in_val, out_val);
                            }
                        }
                    }

                    for other_id in inner_ids.iter().skip(num) {
                        if other_id.is_none() {
                            break;
                        }
                        let (in_val, out_val) = (in_id.id().unwrap(), other_id.id().unwrap());
                        if in_val != out_val {
                            Self::overlap(arc_pcls.as_ref(), in_val, out_val);
                        }
                    }
                }
            }
        }
    }

    #[inline(never)]
    pub fn get_drawable(&self) -> impl Iterator<Item = (f32, f32, f32)> + '_ {
        self.pcls.get_drawable()
    }

    #[inline(never)]
    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, mass: f32, charge: f32) {
        let index = Arc::get_mut(&mut self.pcls)
            .unwrap()
            .push((x, y, radius, mass, charge));
        self.grid.try_insert(index, x, y);
    }

    #[inline(never)]
    pub fn clear(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().clear();
        self.grid.map.clear();
    }
    #[inline(never)]
    pub fn toggle_gravity(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().g_toward_center = !self.pcls.g_toward_center;
    }

    #[inline(never)]
    fn apply_gravity(p: &Particles) {
        for i in 0..p.count {
            if p.g_toward_center {
                let x = p.get_x(i);
                let y = p.get_y(i);
                let r2 = x.abs().powi(2) + y.abs().powi(2);
                let v_x = x / (r2.sqrt());
                let v_y = y / (r2.sqrt());
                p.set_ax(i, -GRAVITY * v_x * (1.0 / (r2 + ANTI_BLACK_HOLE)));
                p.set_ay(i, -GRAVITY * v_y * (1.0 / (r2 + ANTI_BLACK_HOLE)));
            } else {
                p.set_ay(i, -GRAVITY);
            }
        }
    }

    #[inline(never)]
    pub fn apply_coulomb(grid: &Arc<Grid>, pcls: &Arc<Particles>, n_threads: usize) {}

    #[inline(never)]
    fn overlap(p: &Particles, i: usize, j: usize) {
        let xi = p.get_x(i);
        let xj = p.get_x(j);
        let yi = p.get_y(i);
        let yj = p.get_y(j);
        let dx = xi - xj;
        let dy = yi - yj;
        let distance_sq = dx * dx + dy * dy;
        let ri = p.get_r(i);
        let rj = p.get_r(j);
        if distance_sq > (ri + rj) * (ri + rj) {
            return;
        }

        let distance = distance_sq.sqrt();
        let overlap_distance = ri + rj - distance;
        let (normal_x, normal_y) = if distance != 0.0 {
            (dx / distance, dy / distance)
        } else {
            let theta = rand::thread_rng().gen_range(0.0..std::f32::consts::PI * 2.0);
            (theta.cos(), theta.sin())
        };

        let mi = p.get_m(i);
        let mj = p.get_m(j);
        let mass_ratio_1 = mi / (mi + mj);
        let mass_ratio_2 = mj / (mi + mj);

        let correction = RESTITUTION * overlap_distance * 0.5;
        let correction_x = correction * normal_x;
        let correction_y = correction * normal_y;

        p.set_x(i, xi + correction_x * mass_ratio_2);
        p.set_y(i, yi + correction_y * mass_ratio_2);
        p.set_x(j, xj - correction_x * mass_ratio_1);
        p.set_y(j, yj - correction_y * mass_ratio_1);
    }

    #[inline(never)]
    pub fn verlet(p: &Particles, dt: f32) {
        for i in 0..p.count {
            let x = p.get_x(i);
            let y = p.get_y(i);
            let vx = (x - p.get_ox(i)).clamp(-MAX_V, MAX_V) * VELOCITY_DAMPING;
            let vy = (y - p.get_oy(i)).clamp(-MAX_V, MAX_V) * VELOCITY_DAMPING;
            p.set_ox(i, x);
            p.set_oy(i, y);
            p.set_x(i, x + vx + p.get_ax(i) * dt * dt);
            p.set_y(i, y + vy + p.get_ay(i) * dt * dt);
            p.set_ax(i, 0.0);
            p.set_ay(i, 0.0);
        }
    }

    #[inline(never)]
    pub fn constrain(p: &Particles) {
        for i in 0..p.count {
            let x = p.get_x(i);
            let y = p.get_y(i);
            let r = p.get_r(i);
            if WASHING_MACHINE {
                let center_dist = (x * x + y * y).sqrt();
                let factor = if center_dist + r > 1.0 {
                    (1.0 - r) / center_dist
                } else if center_dist - r < 0.3 {
                    (0.3 + r) / center_dist
                } else {
                    1.0
                };
                p.set_x(i, x * factor);
                p.set_y(i, y * factor);
            } else {
                if x + r > 1.0 {
                    p.set_x(i, 1.0 - r);
                } else if x - r < -1.0 {
                    p.set_x(i, -1.0 + r);
                }
                if y + r > 1.0 {
                    p.set_y(i, 1.0 - r);
                } else if y - r < -1.0 {
                    p.set_y(i, -1.0 + r);
                }
            }
        }
    }

    #[inline(never)]
    pub fn stop(&mut self) {
        for i in 0..self.pcls.count {
            self.pcls.set_ox(i, self.pcls.get_x(i));
            self.pcls.set_oy(i, self.pcls.get_y(i));
        }
    }
}
