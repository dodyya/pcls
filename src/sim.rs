use crate::grid::Grid;
use crate::maybe_id::MaybeID;
use crate::particles::Particles;
use crate::types::{Real, TAU};
use rand::Rng;
use rayon::prelude::*;
use std::sync::Arc;

const DT: Real = 1.0 / 60.0;
const SUBSTEPS: usize = 12; // Substeps per step() call.
const GRAVITY: Real = 1.0;

const ANTI_BHOLE: Real = 0.5; // Avoid black holes at center in donut mode
const RESTITUTION: Real = 1.0; // How hard particles bounce off each other, 0.0-1.0
const MAX_V: Real = 0.01; // Maximum velocity restriction
const GRID_DEPTH: usize = 3; // Max. particles per grid cell to process. 3 is reasonable lower limit
const VELOCITY_DAMPING: Real = 0.999999; // Velocity damping per Verlet step
const K: Real = 0.002; // Hue-force strength
const HUE_FORCE_RADIUS: i32 = 3; // Grid "radius" for the hue force
const HOLE_SIZE: Real = 0.3;
const DONUT_SIZE: Real = 1.0;

#[derive(Debug)]
pub struct Simulation {
    pub pcls: Arc<Particles>,
    pub grid: Arc<Grid>,
}

impl Simulation {
    pub fn new(cell_size: Real) -> Self {
        let grid = Grid::new(cell_size, GRID_DEPTH);
        Self {
            pcls: Arc::new(Particles::new(1)),
            grid: Arc::new(grid),
        }
    }
    pub fn step(&mut self) {
        for _ in 0..SUBSTEPS {
            Self::apply_gravity(&self.pcls);
            if self.pcls.hue_force_enabled {
                Self::apply_hue_force(&self.grid, &self.pcls);
            }
            Self::constrain(&self.pcls);
            self.grid.update(&self.pcls);
            Self::resolve_overlaps(&self.grid, &self.pcls);
            Self::verlet(&self.pcls, DT / SUBSTEPS as Real);
        }
    }

    pub fn resolve_overlaps(grid: &Grid, pcls: &Particles) {
        let c = grid.cell_count;
        // Two passes (even/odd columns) so parallel tasks don't write to adjacent columns
        (0..(c + 1) / 2).into_par_iter().for_each(|half_i| {
            Self::overlap_column(half_i * 2, c, grid, pcls);
        });
        (0..c / 2).into_par_iter().for_each(|half_i| {
            Self::overlap_column(half_i * 2 + 1, c, grid, pcls);
        });
    }

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

    fn hue_force_neighbors<'a>(
        out: &mut Vec<&'a [MaybeID]>,
        grid: &'a Grid,
        i: usize,
        j: usize,
        c: usize,
    ) {
        out.clear();
        // Forward-only half-plane (same pattern as `neighbors`, scaled to radius)
        // so each cross-cell pair is visited exactly once.
        for di in -HUE_FORCE_RADIUS..=HUE_FORCE_RADIUS {
            for dj in -HUE_FORCE_RADIUS..=HUE_FORCE_RADIUS {
                if !(dj > 0 || (dj == 0 && di > 0)) {
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

    fn overlap_column(i: usize, c: usize, grid: &Grid, pcls: &Particles) {
        let mut neigh: Vec<&[MaybeID]> = Vec::with_capacity(4);

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

    fn hue_force_column(i: usize, c: usize, grid: &Grid, pcls: &Particles) {
        let mut neigh: Vec<&[MaybeID]> = Vec::with_capacity(24);

        for j in (0..c).rev() {
            let inner = &grid.map[(i, j)];
            Self::hue_force_neighbors(&mut neigh, grid, i, j, c);

            for (idx, in_id) in inner.iter().take_while(|x| x.is_some()).enumerate() {
                let in_val = in_id.unchecked_id();

                for ids in &neigh {
                    for out_id in (*ids).iter().take_while(|x| x.is_some()) {
                        Self::hue_force(pcls, in_val, out_id.unchecked_id());
                    }
                }

                for other_id in inner.iter().skip(idx + 1).take_while(|x| x.is_some()) {
                    Self::hue_force(pcls, in_val, other_id.unchecked_id());
                }
            }
        }
    }

    pub fn get_drawable(&self) -> impl Iterator<Item = (Real, Real, Real, Real)> + '_ {
        self.pcls.get_drawable()
    }

    pub fn is_hue_force_enabled(&self) -> bool {
        self.pcls.hue_force_enabled
    }

    pub fn add_particle(&mut self, x: Real, y: Real, radius: Real, mass: Real, hue: Real) {
        let index = Arc::get_mut(&mut self.pcls)
            .unwrap()
            .push((x, y, radius, mass, hue));
        self.grid.try_insert(index, x, y);
    }

    pub fn clear(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().clear();
        self.grid.map.clear();
    }

    pub fn toggle_gravity(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().g_toward_center = !self.pcls.g_toward_center;
    }

    pub fn toggle_hue_force(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().hue_force_enabled = !self.pcls.hue_force_enabled;
    }

    pub fn toggle_donut(&mut self) {
        Arc::get_mut(&mut self.pcls).unwrap().donut_enabled = !self.pcls.donut_enabled;
    }

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

    pub fn apply_hue_force(grid: &Grid, pcls: &Particles) {
        let c = grid.cell_count;
        (0..c).into_par_iter().for_each(|i| {
            Self::hue_force_column(i, c, grid, pcls);
        });
    }

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
            let theta = rand::thread_rng().gen_range(0.0..TAU);
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

    // Hue-force: like-hues attract, opposite hues (Δh = 0.5) repel; neutral at Δh = 0.25.
    // `-cos(2π · Δh)` is naturally periodic in hue, so wrap-around is free.
    fn hue_force(p: &Particles, i: usize, j: usize) {
        let xi = p.get_x(i);
        let yi = p.get_y(i);
        let hi = p.get_hue(i);
        let mi = p.get_m(i);
        let xj = p.get_x(j);
        let yj = p.get_y(j);
        let hj = p.get_hue(j);
        let mj = p.get_m(j);
        let dx = xi - xj;
        let dy = yi - yj;
        let distance_sq = dx * dx + dy * dy;
        let hue_dist = (hi - hj).abs().min((hi - hj - 1.).abs()).min((hi - hj + 1.).abs());
        // Linear: force ∝ (0.25 - hue_dist), negative (attract) for like-hues,
        // positive (repel) for opposites, zero at quarter-turn. Amplitude is 0.25
        // (vs 1.0 for the cos version) so K is effectively ~4× weaker here.
        // Scale by (touching-distance)² so the force at touching distance is K * factor,
        // independent of particle size.
        let touching = 2.0 * p.min_radius;
        let force = K * (hue_dist - 0.25) * touching * touching / distance_sq;
        let distance = distance_sq.sqrt().max(2.*p.min_radius);

        let fx = force * dx / distance;
        let fy = force * dy / distance;

        // fetch_add (not load+store) — parallel column workers can hit the same
        // particle since the 3-cell reach overlaps adjacent par-iter tasks.
        p.add_ax(i, fx / mi);
        p.add_ay(i, fy / mi);
        p.add_ax(j, -fx / mj);
        p.add_ay(j, -fy / mj);
    }

    pub fn verlet(p: &Particles, dt: Real) {
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

    pub fn constrain(p: &Particles) {
        for i in 0..p.count {
            let x = p.get_x(i);
            let y = p.get_y(i);
            let r = p.get_r(i);
            if p.donut_enabled {
                let center_dist = (x * x + y * y).sqrt();
                let factor = if center_dist + r > DONUT_SIZE {
                    (DONUT_SIZE - r) / center_dist
                } else if center_dist - r < HOLE_SIZE {
                    (HOLE_SIZE + r) / center_dist
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

    pub fn stop(&mut self) {
        for i in 0..self.pcls.count {
            self.pcls.set_ox(i, self.pcls.get_x(i));
            self.pcls.set_oy(i, self.pcls.get_y(i));
        }
    }
}
