use crate::grid::Grid;
use crate::maybe_id::MaybeID;
use crate::particles::Particles;
use rand::Rng;
use std::ops::Range;
use std::sync::Arc;
use std::thread;

const DT: f32 = 1.0 / 60.0;
const SUBSTEPS: usize = 12; // Substeps per step() call.
const NUM_THREADS: usize = 8; // Threads used in resolve_overlaps
const GRAVITY: f32 = 1.0;

const ANTI_BHOLE: f32 = 0.5; // Avoid black holes at center in donut mode
const RESTITUTION: f32 = 1.0; // How hard particles bounce off each other, 0.0-1.0
const MAX_V: f32 = 0.01; // Maximum velocity restriction
const GRID_DEPTH: usize = 3; // Max. particles per grid cell to process. 3 is reasonable lower limit
const VELOCITY_DAMPING: f32 = 0.999999; // Velocity damping per Verlet step
const K: f32 = 0.000000025; //Coulomb's constant
const COULOMB_RADIUS: i32 = 3; //Grid "radius" for Coulomb
const EPSILON: f32 = 0.01; //Coulomb minimum distance

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
    #[inline(always)]
    pub fn step(&mut self) {
        for _ in 0..SUBSTEPS {
            Self::apply_gravity(&self.pcls);
            if self.pcls.coulomb_enabled {
                Self::apply_coulomb(&self.grid, &self.pcls, NUM_THREADS);
            }
            Self::constrain(&self.pcls);
            self.grid.update(&self.pcls);
            Self::resolve_overlaps(&self.grid, &self.pcls, NUM_THREADS);
            Self::verlet(&self.pcls, DT / SUBSTEPS as f32);
        }
    }

    #[inline(always)]
    pub fn resolve_overlaps(grid: &Arc<Grid>, pcls: &Arc<Particles>, n_threads: usize) {
        let c = grid.cell_count;
        if c % n_threads != 0 {
            panic!("Cells to a side is not divisible by number of threads");
        }

        let thread_width = c / n_threads;
        thread::scope(|s| {
            for n in 0..n_threads {
                let range = n * thread_width..n * thread_width + thread_width / 2;
                s.spawn(move || Self::overlap_chunk(range, c, Arc::clone(grid), Arc::clone(pcls)));
            }
        });

        thread::scope(|s| {
            for n in 0..n_threads {
                let range = n * thread_width + thread_width / 2..(n + 1) * thread_width;
                s.spawn(move || Self::overlap_chunk(range, c, Arc::clone(grid), Arc::clone(pcls)));
            }
        });
    }

    #[inline(always)]
    fn neighbors<'a>(out: &mut Vec<&'a [MaybeID]>, grid: &'a Grid, i: usize, j: usize, c: usize) {
        out.clear();
        if i > 0 && j + 1 < c {
            out.push(&grid.map[(i - 1, j + 1)]);
        }
        if i + 1 < c {
            out.push(&grid.map[(i + 1, j)]);
        }
        if j + 1 < c {
            out.push(&grid.map[(i, j + 1)]);
        }
        if i + 1 < c && j + 1 < c {
            out.push(&grid.map[(i + 1, j + 1)]);
        }
    }

    #[inline(always)]
    fn coulomb_neighbors<'a>(
        out: &mut Vec<&'a [MaybeID]>,
        grid: &'a Grid,
        i: usize,
        j: usize,
        c: usize,
    ) {
        out.clear();
        for di in -COULOMB_RADIUS..=COULOMB_RADIUS {
            for dj in -COULOMB_RADIUS..=COULOMB_RADIUS {
                if di == 0 && dj == 0 {
                    continue;
                }
                let ni = i as i32 + di;
                let nj = j as i32 + dj;
                if ni >= 0 && ni < c as i32 && nj >= 0 && nj < c as i32 {
                    out.push(&grid.map[(ni as usize, nj as usize)]);
                }
            }
        }
    }

    #[inline(always)]
    fn overlap_chunk(
        x_range: Range<usize>,
        c: usize,
        arc_grid: Arc<Grid>,
        arc_pcls: Arc<Particles>,
    ) {
        let grid = arc_grid.as_ref();
        let pcls = arc_pcls.as_ref();

        let mut neigh: Vec<&[MaybeID]> = Vec::with_capacity(4);

        for i in x_range {
            for j in (0..c).rev() {
                let inner = &grid.map[(i, j)];
                Self::neighbors(&mut neigh, grid, i, j, c);

                for (idx, in_id) in inner.iter().take_while(|x| x.is_some()).enumerate() {
                    let in_val = in_id.unchecked_id();

                    for ids in &neigh {
                        for out_id in (*ids).iter().take_while(|x| x.is_some()) {
                            Self::overlap(pcls, in_val, out_id.unchecked_id());
                        }
                    }

                    for other_id in inner.iter().skip(idx + 1).take_while(|x| x.is_some()) {
                        Self::overlap(pcls, in_val, other_id.unchecked_id());
                    }
                }
            }
        }
    }

    #[inline(always)]
    fn coulomb_chunk(
        x_range: Range<usize>,
        c: usize,
        arc_grid: Arc<Grid>,
        arc_pcls: Arc<Particles>,
    ) {
        let grid = arc_grid.as_ref();
        let pcls = arc_pcls.as_ref();

        let mut magnetic_neigh: Vec<&[MaybeID]> = Vec::with_capacity(24);

        for i in x_range {
            for j in (0..c).rev() {
                let inner = &grid.map[(i, j)];
                Self::coulomb_neighbors(&mut magnetic_neigh, grid, i, j, c);

                for (idx, in_id) in inner.iter().take_while(|x| x.is_some()).enumerate() {
                    let in_val = in_id.unchecked_id();

                    for ids in &magnetic_neigh {
                        for out_id in (*ids).iter().take_while(|x| x.is_some()) {
                            Self::coulomb(pcls, in_val, out_id.unchecked_id());
                        }
                    }

                    for other_id in inner.iter().skip(idx + 1).take_while(|x| x.is_some()) {
                        Self::coulomb(pcls, in_val, other_id.unchecked_id());
                    }
                }
            }
        }
    }

    #[inline(always)]
    pub fn get_drawable(&self) -> impl Iterator<Item = (f32, f32, f32, f32)> + '_ {
        self.pcls.get_drawable()
    }

    #[inline(always)]
    pub fn is_coulomb_enabled(&self) -> bool {
        self.pcls.coulomb_enabled
    }

    #[inline(always)]
    pub fn add_particle(&mut self, x: f32, y: f32, radius: f32, mass: f32, charge: f32) {
        let index = Arc::get_mut(&mut self.pcls)
            .unwrap()
            .push((x, y, radius, mass, charge));
        self.grid.try_insert(index, x, y);
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().clear();
        self.grid.map.clear();
    }

    #[inline(always)]
    pub fn toggle_gravity(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().g_toward_center = !self.pcls.g_toward_center;
    }

    #[inline(always)]
    pub fn toggle_coulomb(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().coulomb_enabled = !self.pcls.coulomb_enabled;
    }

    #[inline(always)]
    pub fn toggle_donut(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().donut_enabled = !self.pcls.donut_enabled;
    }

    #[inline(always)]
    fn apply_gravity(p: &Particles) {
        for i in 0..p.count {
            if p.g_toward_center {
                let x = p.get_x(i);
                let y = p.get_y(i);
                let r2 = x.abs().powi(2) + y.abs().powi(2);
                let v_x = x / (r2.sqrt());
                let v_y = y / (r2.sqrt());
                p.set_ax(i, -GRAVITY * v_x * (1.0 / (r2 + ANTI_BHOLE)));
                p.set_ay(i, -GRAVITY * v_y * (1.0 / (r2 + ANTI_BHOLE)));
            } else {
                p.set_ay(i, -GRAVITY);
            }
        }
    }

    #[inline(always)]
    pub fn apply_coulomb(grid: &Arc<Grid>, pcls: &Arc<Particles>, n_threads: usize) {
        let c = grid.cell_count;
        if c % n_threads != 0 {
            panic!("Cells to a side is not divisible by number of threads");
        }

        let thread_width = c / n_threads;
        thread::scope(|s| {
            for n in 0..n_threads {
                let range = n * thread_width..n * thread_width + thread_width;
                s.spawn(move || Self::coulomb_chunk(range, c, Arc::clone(grid), Arc::clone(pcls)));
            }
        });
    }

    #[inline(always)]
    fn overlap(p: &Particles, i: usize, j: usize) {
        let xi = p.get_x(i);
        let yi = p.get_y(i);
        let ri = p.get_r(i);
        let xj = p.get_x(j);
        let yj = p.get_y(j);
        let rj = p.get_r(j);
        let dx = xi - xj;
        let dy = yi - yj;
        let distance_sq = dx * dx + dy * dy;
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

    fn coulomb(p: &Particles, i: usize, j: usize) {
        let xi = p.get_x(i);
        let yi = p.get_y(i);
        let ci = p.get_c(i);
        let mi = p.get_m(i);
        let xj = p.get_x(j);
        let yj = p.get_y(j);
        let cj = p.get_c(j);
        let mj = p.get_m(j);
        let dx = xi - xj;
        let dy = yi - yj;
        let distance_sq = dx * dx + dy * dy;
        let force = K * ci * cj / (distance_sq);
        let distance = distance_sq.sqrt().max(EPSILON);

        let fx = force * dx / distance;
        let fy = force * dy / distance;

        p.set_ax(i, p.get_ax(i) + fx / mi);
        p.set_ay(i, p.get_ay(i) + fy / mi);
        p.set_ax(j, p.get_ax(j) - fx / mj);
        p.set_ay(j, p.get_ay(j) - fy / mj);
    }

    #[inline(always)]
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

    #[inline(always)]
    pub fn constrain(p: &Particles) {
        for i in 0..p.count {
            let x = p.get_x(i);
            let y = p.get_y(i);
            let r = p.get_r(i);
            if p.donut_enabled {
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

    #[inline(always)]
    pub fn stop(&mut self) {
        for i in 0..self.pcls.count {
            self.pcls.set_ox(i, self.pcls.get_x(i));
            self.pcls.set_oy(i, self.pcls.get_y(i));
        }
    }
}
