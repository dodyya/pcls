use crate::array::Array3D;

use crate::maybe_id::MaybeID;
use crate::particles::{Particles, CELL_UNREGISTERED};
use crate::types::Real;

#[derive(Debug)]
pub struct Grid {
    pub cell_count: usize,
    pub depth: usize,
    pub map: Array3D<MaybeID>,
}

impl Grid {
    pub fn new(cell_size: Real, depth: usize) -> Self {
        let cell_count = (2.0 / cell_size).ceil() as usize;
        Self {
            cell_count,
            depth,
            map: Array3D::default(cell_count, cell_count, depth),
        }
    }

    pub fn index(&self, x: Real, y: Real) -> (usize, usize) {
        let cell_size = 2.0 / self.cell_count as Real;
        let x_index = ((x + 1.0) / cell_size).floor() as usize;
        let y_index = ((y + 1.0) / cell_size).floor() as usize;
        (x_index, y_index)
    }

    // Insert a freshly-pushed particle. Used by add_particle; safe to call once
    // per particle right after Particles::push.
    pub fn try_insert(&self, particles: &Particles, id: usize, x: Real, y: Real) {
        let ind = self.index(x, y);
        if self.in_bounds(ind) && self.ind_insert(id, ind) {
            particles.set_cell(id, ind);
        }
        // else: stays CELL_UNREGISTERED — the next update() will retry.
    }

    // Returns true if a free slot was found in the cell.
    pub fn ind_insert(&self, id: usize, ind: (usize, usize)) -> bool {
        let cell = &self.map[ind];
        for j in 0..self.depth {
            if cell[j].is_none() {
                cell[j].set(id);
                return true;
            }
        }
        false
    }

    // Removes id while keeping the cell densely packed: empty slots only at
    // the tail. The overlap/hue-force loops in sim.rs walk cells with
    // `take_while(is_some)` and would skip past anything behind a hole, so
    // leaving a None in the middle hides every particle behind it.
    pub fn ind_remove(&self, id: usize, ind: (usize, usize)) {
        let cell = &self.map[ind];
        let mut found: Option<usize> = None;
        let mut last: Option<usize> = None;
        for j in 0..self.depth {
            if cell[j].is_some() {
                if cell[j].unchecked_id() == id {
                    found = Some(j);
                }
                last = Some(j);
            } else {
                break;
            }
        }
        let (Some(j), Some(k)) = (found, last) else {
            // Invariant: if a particle's cell_xs/cell_ys says it's here, it is.
            // Tripping this means cell tracking and the grid drifted apart.
            debug_assert!(false, "ind_remove: id {} not found in cell {:?}", id, ind);
            return;
        };
        if j == k {
            cell[j].make_empty();
        } else {
            cell[j].set(cell[k].unchecked_id());
            cell[k].make_empty();
        }
    }

    pub fn in_bounds(&self, (ix, iy): (usize, usize)) -> bool {
        ix < self.cell_count && iy < self.cell_count
    }

    // Incremental update: for each particle, compare its current (x, y)
    // against the cell it's registered in. Migrate only when the cell changes.
    // Safe to call repeatedly; relies on the per-particle cell field being the
    // single source of truth for "where am I in the grid".
    pub fn update(&self, particles: &Particles) {
        for i in 0..particles.count {
            let x = particles.get_x(i);
            let y = particles.get_y(i);
            let new = self.index(x, y);
            if !self.in_bounds(new) {
                continue;
            }

            let (ox, oy) = particles.get_cell(i);
            let registered = ox != CELL_UNREGISTERED;

            if registered && (ox, oy) == new {
                continue;
            }

            // Remove from the old cell first, then attempt the new one. On
            // failure mark the particle unregistered: gravity + verlet keep
            // moving it, and next substep we retry insertion at its (possibly
            // different) current cell. This matches the rebuild-from-scratch
            // semantics where overflow particles silently sit out a substep
            // rather than getting pinned to a stale neighborhood — the latter
            // was self-reinforcing in hue-attractive dense clusters.
            if registered {
                self.ind_remove(i, (ox, oy));
            }
            if self.ind_insert(i, new) {
                particles.set_cell(i, new);
            } else {
                particles.set_cell_unregistered(i);
            }
        }
    }
}
