use crate::grid::Grid;
use crate::maybe_id::MaybeID;
use crate::particles::Particles;
use crate::types::{Real, TAU};
use rand::Rng;
use rayon::prelude::*;
use std::sync::Arc;

pub const PARTICLE_RADIUS: Real = 1.0 / 300.0;
const PARTICLE_MASS: Real = 1.0;

const DT: Real = 1.0 / 60.0;
const SUBSTEPS: usize = 12; // Substeps per step() call.
const GRAVITY: Real = 1.0;

const ANTI_BHOLE: Real = 0.5; // Avoid black holes at center in donut mode
const RESTITUTION: Real = 0.5; // How hard particles bounce off each other, 0.0-1.0
const MAX_V: Real = 0.0025; // Maximum velocity restriction
const GRID_DEPTH: usize = 4; // Max. particles per grid cell to process. 3 is reasonable lower limit
const VELOCITY_DAMPING: Real = 0.999999; // Velocity damping per Verlet step
const K: Real = 300.*256.*PARTICLE_RADIUS;//*(1./256./PARTICLE_RADIUS); // Hue-force strength
const HUE_FORCE_RADIUS: i32 = 3; // Grid "radius" for the hue force
const HOLE_SIZE: Real = 0.05;
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
            pcls: Arc::new(Particles::new(60_000)),
            grid: Arc::new(grid),
        }
    }
    pub fn step(&mut self) {
        // BENCH: per-phase timing, accumulated across all substeps, printed every 60 frames.
        // Remove this block to restore the plain loop.
        use std::sync::atomic::{AtomicU64, Ordering::Relaxed as R};
        use std::time::Instant;
        static FRAMES: AtomicU64 = AtomicU64::new(0);
        static T_GRAV: AtomicU64 = AtomicU64::new(0);
        static T_HUE: AtomicU64 = AtomicU64::new(0);
        static T_CON: AtomicU64 = AtomicU64::new(0);
        static T_GRID: AtomicU64 = AtomicU64::new(0);
        static T_OVL: AtomicU64 = AtomicU64::new(0);
        static T_VER: AtomicU64 = AtomicU64::new(0);

        for i in 0..SUBSTEPS {
            let t0 = Instant::now();
            Self::apply_gravity(&self.pcls);
            let t1 = Instant::now();
            if self.pcls.hue_force_enabled && i%4==0 {
                Self::apply_hue_force(&self.grid, &self.pcls);
            }
            let t2 = Instant::now();
            Self::constrain(&self.pcls);
            let t3 = Instant::now();
            self.grid.update(&self.pcls);
            let t4 = Instant::now();
            Self::resolve_overlaps(&self.grid, &self.pcls);
            let t5 = Instant::now();
            Self::verlet(&self.pcls, DT / SUBSTEPS as Real);
            let t6 = Instant::now();
            T_GRAV.fetch_add((t1 - t0).as_nanos() as u64, R);
            T_HUE.fetch_add((t2 - t1).as_nanos() as u64, R);
            T_CON.fetch_add((t3 - t2).as_nanos() as u64, R);
            T_GRID.fetch_add((t4 - t3).as_nanos() as u64, R);
            T_OVL.fetch_add((t5 - t4).as_nanos() as u64, R);
            T_VER.fetch_add((t6 - t5).as_nanos() as u64, R);
        }

        let n = FRAMES.fetch_add(1, R) + 1;
        if n % 60 == 0 {
            // Wall-clock FPS over the last 60-frame window (includes render+swap,
            // not just the sim phases below) so framerate drops show up here.
            static LAST_PRINT: std::sync::Mutex<Option<Instant>> = std::sync::Mutex::new(None);
            let now = Instant::now();
            let fps = {
                let mut lp = LAST_PRINT.lock().unwrap();
                let fps = lp.map_or(0.0, |prev| 60.0 / (now - prev).as_secs_f64());
                *lp = Some(now);
                fps
            };
            let f = (n * SUBSTEPS as u64) as f64 * 1000.0; // ns→µs / substep avg
            eprintln!(
                "[sim n={:>5} pcls f={:>4}] fps={:>5.1} grav={:>5.1} hue={:>6.1} con={:>5.1} grid={:>5.1} ovl={:>6.1} ver={:>5.1}   (µs/substep avg, hue on={})",
                self.pcls.count,
                n,
                fps,
                T_GRAV.load(R) as f64 / f,
                T_HUE.load(R) as f64 / f,
                T_CON.load(R) as f64 / f,
                T_GRID.load(R) as f64 / f,
                T_OVL.load(R) as f64 / f,
                T_VER.load(R) as f64 / f,
                self.pcls.hue_force_enabled,
            );
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

    fn hue_force_column(
        i: usize,
        c: usize,
        grid: &Grid,
        pcls: &Particles,
        accum: &mut [(Real, Real)],
    ) {
        let mut neigh: Vec<&[MaybeID]> = Vec::with_capacity(24);

        for j in (0..c).rev() {
            let inner = &grid.map[(i, j)];
            Self::hue_force_neighbors(&mut neigh, grid, i, j, c);

            for (idx, in_id) in inner.iter().take_while(|x| x.is_some()).enumerate() {
                let in_val = in_id.unchecked_id();

                for ids in &neigh {
                    for out_id in (*ids).iter().take_while(|x| x.is_some()) {
                        Self::hue_force(pcls, accum, in_val, out_id.unchecked_id());
                    }
                }

                for other_id in inner.iter().skip(idx + 1).take_while(|x| x.is_some()) {
                    Self::hue_force(pcls, accum, in_val, other_id.unchecked_id());
                }
            }
        }
    }

    pub fn get_drawable(&self) -> impl Iterator<Item = (Real, Real, Real)> + '_ {
        self.pcls.get_drawable()
    }

    pub fn is_hue_force_enabled(&self) -> bool {
        self.pcls.hue_force_enabled
    }

    pub fn add_particle(&mut self, x: Real, y: Real, hue: Real) {
        let index = Arc::get_mut(&mut self.pcls).unwrap().push((x, y, hue));
        self.grid.try_insert(&self.pcls, index, x, y);
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
        let n = pcls.count;

        // Each rayon worker folds into a thread-local plain accumulator
        // (no atomics on the hot path); the per-thread buffers are then
        // reduced and applied to the canonical axs/ays once.
        let totals = (0..c)
            .into_par_iter()
            .fold(
                || vec![(0.0 as Real, 0.0 as Real); n],
                |mut local, col| {
                    Self::hue_force_column(col, c, grid, pcls, &mut local);
                    local
                },
            )
            .reduce(
                || vec![(0.0 as Real, 0.0 as Real); n],
                |mut a, b| {
                    for k in 0..n {
                        a[k].0 += b[k].0;
                        a[k].1 += b[k].1;
                    }
                    a
                },
            );

        for k in 0..n {
            pcls.add_ax(k, totals[k].0);
            pcls.add_ay(k, totals[k].1);
        }
    }

    fn overlap(p: &Particles, i: usize, j: usize) {
        let xi = p.get_x(i);
        let yi = p.get_y(i);
        let xj = p.get_x(j);
        let yj = p.get_y(j);
        let dx = xi - xj;
        let dy = yi - yj;
        let contact = 2.0 * PARTICLE_RADIUS;
        let distance_sq = dx * dx + dy * dy;
        if distance_sq > contact * contact {
            return;
        }

        let distance = distance_sq.sqrt();
        let overlap_distance = contact - distance;
        let (normal_x, normal_y) = if distance != 0.0 {
            (dx / distance, dy / distance)
        } else {
            let theta = rand::thread_rng().gen_range(0.0..TAU);
            (theta.cos(), theta.sin())
        };

        // Equal masses → each particle takes half the correction.
        let correction = RESTITUTION * overlap_distance * 0.25;
        let correction_x = correction * normal_x;
        let correction_y = correction * normal_y;

        p.set_x(i, xi + correction_x);
        p.set_y(i, yi + correction_y);
        p.set_x(j, xj - correction_x);
        p.set_y(j, yj - correction_y);
    }

    // Hue-force: like-hues attract, opposite hues (Δh = 0.5) repel; neutral at Δh = 0.25.
    // `-cos(2π · Δh)` is naturally periodic in hue, so wrap-around is free.
    fn hue_force(p: &Particles, accum: &mut [(Real, Real)], i: usize, j: usize) {
        let xi = p.get_x(i);
        let yi = p.get_y(i);
        let hi = p.get_hue(i);
        let xj = p.get_x(j);
        let yj = p.get_y(j);
        let hj = p.get_hue(j);
        let dx = xi - xj;
        let dy = yi - yj;

        let hue_dist = (hi - hj)
            .abs()
            .min((hi - hj - 1.).abs())
            .min((hi - hj + 1.).abs());
        // Linear: force ∝ (0.25 - hue_dist), negative (attract) for like-hues,
        // positive (repel) for opposites, zero at quarter-turn. Amplitude is 0.25
        // (vs 1.0 for the cos version) so K is effectively ~4× weaker here.
        // Floor the distance at the contact distance: the 1/d² term is capped
        // at the surface so force is K * factor at contact and can't keep
        // growing as particles overlap.
        let contact = 2.0 * PARTICLE_RADIUS;
        let d2 = (dx * dx + dy * dy).max(contact * contact);
        let force = K * (hue_dist - 0.25) * contact * contact / d2;
        let distance = d2.sqrt();

        let fx = force * dx / distance;
        let fy = force * dy / distance;

        // Plain += into a thread-local accumulator — apply_hue_force reduces
        // the per-thread buffers and folds them into the canonical axs/ays.
        accum[i].0 += fx / PARTICLE_MASS;
        accum[i].1 += fy / PARTICLE_MASS;
        accum[j].0 -= fx / PARTICLE_MASS;
        accum[j].1 -= fy / PARTICLE_MASS;
    }

    pub fn verlet(p: &Particles, dt: Real) {
        for i in 0..p.count {
            let x = p.get_x(i);
            let y = p.get_y(i);
            let vx = (x - p.get_ox(i)) * VELOCITY_DAMPING;
            let vy = (y - p.get_oy(i)) * VELOCITY_DAMPING;
            p.set_ox(i, x);
            p.set_oy(i, y);
            p.set_x(i, x + (vx + p.get_ax(i) * dt * dt).clamp(-MAX_V, MAX_V));
            p.set_y(i, y + (vy + p.get_ay(i) * dt * dt).clamp(-MAX_V, MAX_V));
            p.set_ax(i, 0.0);
            p.set_ay(i, 0.0);
        }
    }

    pub fn constrain(p: &Particles) {
        let r = PARTICLE_RADIUS;
        for i in 0..p.count {
            let x = p.get_x(i);
            let y = p.get_y(i);
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
