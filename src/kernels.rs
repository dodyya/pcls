// Per substep: grid_build → gravity → (hue_force) → constrain → grid_build → overlap → verlet.
// Positions (XyHue) and integration state (RawParticle) live in two parallel buffers so
// the neighbour-gather kernels (hue_force, overlap, grid_build) read only the 12 B of
// XyHue they actually use instead of 28 B of mixed AoS state, and so hue_force can use
// the blessed `&[T]` + DisjointSlice idiom (read positions, write own ax/ay).
//
// overlap still takes `*mut XyHue` — the racy in-place position read/write across
// neighbours is intentional iterative refinement, which DisjointSlice can't express.
// Each thread writes only its own slot, so the aliasing is benign per the audit.
//
// counters[c] may exceed GRID_DEPTH; cell iteration caps at DEPTH and drops overflow.

use crate::particles::{RawParticle, XyHue};
use crate::sim::{
    ANTI_BHOLE, DIAG_CON_FIRES, DIAG_CON_MAG, DIAG_HUE_MAG, DIAG_OVL_FIRES, DIAG_OVL_MAG,
    DIAG_VRP_DISP, DIAG_VRP_FIRES, DIAG_VRP_VEL, DONUT_SIZE, GRAVITY, GRID_DEPTH, HOLE_SIZE,
    HUE_FORCE_RADIUS, K, MAX_OVERLAP_PUSH, MAX_V, PARTICLE_RADIUS, RESTITUTION, VELOCITY_DAMPING,
};
use crate::types::Real;
use cuda_device::atomic::{AtomicOrdering, DeviceAtomicU32};
use cuda_device::{DisjointSlice, kernel, thread};
use cuda_host::cuda_module;

// Newton 1/sqrt. cuda-oxide can't currently lower libdevice's __nv_rsqrtf on this
// toolchain via the embedded-NVVM path, so we hand-roll one. Quake-III bit-twiddle
// initial guess (≈7 bits), then 3 Newton iterations ≈ 22 bits — enough for sim use.
// Callers replace `coord / sqrt(r²)` with `coord * fast_rsqrt(r²)`, which is the
// canonical pattern that maps to a single MUFU.RSQ on NVIDIA hw.
#[inline(always)]
fn fast_rsqrt(x: Real) -> Real {
    if x <= 0.0 {
        return 0.0;
    }
    let bits: u32 = x.to_bits();
    let guess_bits = 0x5f3759df_u32.wrapping_sub(bits >> 1);
    let mut g = Real::from_bits(guess_bits);
    let half_x = 0.5 * x;
    g = g * (1.5 - half_x * g * g);
    g = g * (1.5 - half_x * g * g);
    g = g * (1.5 - half_x * g * g);
    g
}

// sqrt via x * (1/sqrt(x)) — saves the divide that a Newton sqrt would carry.
#[inline(always)]
fn fast_sqrt(x: Real) -> Real {
    x * fast_rsqrt(x)
}

#[inline(always)]
fn fast_abs(x: Real) -> Real {
    if x < 0.0 { -x } else { x }
}

#[inline(always)]
fn fast_min(a: Real, b: Real) -> Real {
    if a < b { a } else { b }
}

#[inline(always)]
fn fast_max(a: Real, b: Real) -> Real {
    if a > b { a } else { b }
}

// f32 bit pattern is monotonic for non-negative floats, so atomic-max on bits tracks
// max magnitude. (DeviceAtomicF32::fetch_max is unavailable — PTX float atomics support
// only load/store/fetch_add/swap.) Slot layout: see DIAG_* consts in sim.rs.
#[inline(always)]
fn diag_atomic(diag: &[u32], slot: usize) -> &DeviceAtomicU32 {
    // Safety: slot < diag.len() asserted by caller; transparent over u32 layout.
    unsafe { DeviceAtomicU32::from_ptr(diag.as_ptr().add(slot) as *mut u32) }
}

#[inline(always)]
fn diag_max_mag(diag: &[u32], slot: usize, val: Real) {
    let bits = (fast_abs(val) as f32).to_bits();
    diag_atomic(diag, slot).fetch_max(bits, AtomicOrdering::Relaxed);
}

#[inline(always)]
fn diag_incr(diag: &[u32], slot: usize) {
    diag_atomic(diag, slot).fetch_add(1, AtomicOrdering::Relaxed);
}

#[cuda_module]
pub mod kernels {
    use crate::types::Real;

    use super::*;

    // mode: 0 = downward gravity (target ignored); 1 = pull toward (target_x, target_y).
    // The toward-a-point path is shared between space-bar center gravity (target = origin)
    // and shift+left-drag mouse pull (target = cursor in sim space).
    #[kernel]
    pub fn gravity(
        positions: &[XyHue],
        mut integ: DisjointSlice<RawParticle>,
        mode: u32,
        target_x: Real,
        target_y: Real,
    ) {
        let idx = thread::index_1d();
        let i = idx.get();
        if let Some(p) = integ.get_mut(idx) {
            if mode != 0 {
                let pos = positions[i];
                let dx = pos.x - target_x;
                let dy = pos.y - target_y;
                let r2 = dx * dx + dy * dy;
                let inv_r = fast_rsqrt(r2);
                let inv_falloff = 1.0 / (r2 + ANTI_BHOLE);
                p.ax = -GRAVITY * dx * inv_r * inv_falloff;
                p.ay = -GRAVITY * dy * inv_r * inv_falloff;
            } else {
                p.ax = 0.0;
                p.ay = -GRAVITY;
            }
        }
    }

    #[kernel]
    pub fn constrain(mut positions: DisjointSlice<XyHue>, donut_enabled: u32, diag: &[u32]) {
        if let Some(p) = positions.get_mut(thread::index_1d()) {
            let r = PARTICLE_RADIUS;
            let old_x = p.x;
            let old_y = p.y;
            if donut_enabled != 0 {
                let r2 = p.x * p.x + p.y * p.y;
                let inv = fast_rsqrt(r2);
                // sqrt via r2 * (1/sqrt(r2)); inv = 0 at the origin → center_dist = 0.
                let center_dist = r2 * inv;
                let factor = if center_dist + r > DONUT_SIZE {
                    (DONUT_SIZE - r) * inv
                } else if center_dist - r < HOLE_SIZE {
                    (HOLE_SIZE + r) * inv
                } else {
                    1.0
                };
                p.x *= factor;
                p.y *= factor;
            } else {
                if p.x + r > 1.0 {
                    p.x = 1.0 - r;
                } else if p.x - r < -1.0 {
                    p.x = -1.0 + r;
                }
                if p.y + r > 1.0 {
                    p.y = 1.0 - r;
                } else if p.y - r < -1.0 {
                    p.y = -1.0 + r;
                }
            }
            let cdx = p.x - old_x;
            let cdy = p.y - old_y;
            let cmag = fast_sqrt(cdx * cdx + cdy * cdy);
            diag_max_mag(diag, DIAG_CON_MAG, cmag);
            if cmag > MAX_V {
                diag_incr(diag, DIAG_CON_FIRES);
            }
        }
    }

    #[kernel]
    pub fn verlet(
        mut positions: DisjointSlice<XyHue>,
        mut integ: DisjointSlice<RawParticle>,
        dt: Real,
        diag: &[u32],
    ) {
        let pos_opt = positions.get_mut(thread::index_1d());
        let integ_opt = integ.get_mut(thread::index_1d());
        if let (Some(pos), Some(p)) = (pos_opt, integ_opt) {
            let vx = (pos.x - p.ox) * VELOCITY_DAMPING;
            let vy = (pos.y - p.oy) * VELOCITY_DAMPING;
            let carried = fast_sqrt(vx * vx + vy * vy);
            diag_max_mag(diag, DIAG_VRP_VEL, carried);

            p.ox = pos.x;
            p.oy = pos.y;
            let dx_raw = vx + p.ax * dt * dt;
            let dy_raw = vy + p.ay * dt * dt;
            let pre_clamp = fast_sqrt(dx_raw * dx_raw + dy_raw * dy_raw);
            diag_max_mag(diag, DIAG_VRP_DISP, pre_clamp);

            let mut clamped = false;
            let dx = if dx_raw < -MAX_V {
                clamped = true;
                -MAX_V
            } else if dx_raw > MAX_V {
                clamped = true;
                MAX_V
            } else {
                dx_raw
            };
            let dy = if dy_raw < -MAX_V {
                clamped = true;
                -MAX_V
            } else if dy_raw > MAX_V {
                clamped = true;
                MAX_V
            } else {
                dy_raw
            };
            if clamped {
                diag_incr(diag, DIAG_VRP_FIRES);
            }
            pos.x += dx;
            pos.y += dy;
            p.ax = 0.0;
            p.ay = 0.0;
        }
    }

    // Reset encoded velocity (ox = x, oy = y) for every particle.
    #[kernel]
    pub fn zero_velocity(positions: &[XyHue], mut integ: DisjointSlice<RawParticle>) {
        let idx = thread::index_1d();
        let i = idx.get();
        if let Some(p) = integ.get_mut(idx) {
            let pos = positions[i];
            p.ox = pos.x;
            p.oy = pos.y;
        }
    }

    // Caller must zero `counters` before launching. Slots beyond DEPTH are dropped;
    // hue/overlap cap their reads at DEPTH so the overflow is silently ignored.
    //
    // `cells` is `*mut u32` because slot indices come from the atomic counter (dynamic,
    // not thread-witnessed) — DisjointSlice can't model "disjoint by atomic counter".
    #[kernel]
    pub unsafe fn grid_build(
        positions: &[XyHue],
        cells: *mut u32,
        counters: &[u32],
        cell_count: u32,
    ) {
        let i = thread::index_1d().get();
        if i >= positions.len() {
            return;
        }
        let p = positions[i];
        let cell_size = 2.0 / cell_count as Real;
        let cx_f = (p.x + 1.0) / cell_size;
        let cy_f = (p.y + 1.0) / cell_size;
        if cx_f < 0.0
            || cy_f < 0.0
            || cx_f >= cell_count as Real
            || cy_f >= cell_count as Real
        {
            return;
        }
        let cx = cx_f as u32;
        let cy = cy_f as u32;
        let cell_idx = (cx * cell_count + cy) as usize;
        // Safety: cell_idx < counters.len(); atomic ref via interior-mutable transparent type.
        let counter = unsafe {
            DeviceAtomicU32::from_ptr(counters.as_ptr().add(cell_idx) as *mut u32)
        };
        let slot = counter.fetch_add(1, AtomicOrdering::Relaxed);
        if slot < GRID_DEPTH as u32 {
            let slot_idx = cell_idx * GRID_DEPTH + slot as usize;
            // Safety: atomic counter gives each thread a unique slot in [0, GRID_DEPTH),
            // so writes don't race on the same address; `cells` is sized cell_count^2 * GRID_DEPTH.
            unsafe {
                *cells.add(slot_idx) = i as u32;
            }
        }
    }

    // Each thread reads neighbours' (x, y, hue) within HUE_FORCE_RADIUS cells and writes
    // only its own (ax, ay). Positions and integration state are in separate buffers, so
    // the read source and write target are disjoint — this is the blessed safe pattern:
    // `&[T]` read + `DisjointSlice<T>` per-thread write, no `unsafe`.
    #[kernel]
    pub fn hue_force(
        positions: &[XyHue],
        mut integ: DisjointSlice<RawParticle>,
        cells: &[u32],
        counters: &[u32],
        cell_count: u32,
        diag: &[u32],
    ) {
        let idx = thread::index_1d();
        let i = idx.get();
        if i >= positions.len() {
            return;
        }
        let me = positions[i];
        let cell_size = 2.0 / cell_count as Real;
        let mcx = ((me.x + 1.0) / cell_size) as i32;
        let mcy = ((me.y + 1.0) / cell_size) as i32;
        let cc = cell_count as i32;
        let contact = 2.0 * PARTICLE_RADIUS;
        let r = HUE_FORCE_RADIUS;

        let mut fx = 0.0;
        let mut fy = 0.0;

        let mut dy = -r;
        while dy <= r {
            let mut dx = -r;
            while dx <= r {
                let nx = mcx + dx;
                let ny = mcy + dy;
                if !(nx < 0 || nx >= cc || ny < 0 || ny >= cc) {
                    let cell_idx = (nx as u32 * cell_count + ny as u32) as usize;
                    let n_in_cell = counters[cell_idx];
                    let n = if n_in_cell > GRID_DEPTH as u32 {
                        GRID_DEPTH as u32
                    } else {
                        n_in_cell
                    };
                    let mut slot = 0u32;
                    while slot < n {
                        let j = cells[cell_idx * GRID_DEPTH + slot as usize] as usize;
                        if j != i {
                            let other = positions[j];
                            let xdiff = me.x - other.x;
                            let ydiff = me.y - other.y;
                            let hdiff = me.hue - other.hue;
                            let hd0 = fast_abs(hdiff);
                            let hd1 = fast_abs(hdiff - 1.0);
                            let hd2 = fast_abs(hdiff + 1.0);
                            let hue_dist = fast_min(fast_min(hd0, hd1), hd2);
                            let d2_raw = xdiff * xdiff + ydiff * ydiff;
                            let d2 = fast_max(d2_raw, contact * contact);
                            let force = K * (hue_dist - 0.25) * contact * contact / d2;
                            // coord / sqrt(d2) = coord * rsqrt(d2) — one MUFU.RSQ vs sqrt+div.
                            let inv_d = fast_rsqrt(d2);
                            fx += force * xdiff * inv_d;
                            fy += force * ydiff * inv_d;
                        }
                        slot += 1;
                    }
                }
                dx += 1;
            }
            dy += 1;
        }
        let hue_mag = fast_sqrt(fx * fx + fy * fy);
        diag_max_mag(diag, DIAG_HUE_MAG, hue_mag);
        if let Some(p) = integ.get_mut(idx) {
            p.ax += fx;
            p.ay += fy;
        }
    }

    // 3x3-cell overlap pass; pair (i,j) is processed by both threads, each shoving itself
    // by the half-correction → net matches the CPU split-correction update.
    //
    // `positions` is `*mut XyHue` because reads of neighbours' x/y alias writes to own x/y
    // on the same buffer; the racy read is part of the iterative refinement.
    #[kernel]
    pub unsafe fn overlap(
        positions: *mut XyHue,
        cells: &[u32],
        counters: &[u32],
        cell_count: u32,
        count: u32,
        diag: &[u32],
    ) {
        let i = thread::index_1d().get();
        if i >= count as usize {
            return;
        }
        // Safety: i < count and positions is sized to count.
        let me = unsafe { *positions.add(i) };
        let cell_size = 2.0 / cell_count as Real;
        let mcx = ((me.x + 1.0) / cell_size) as i32;
        let mcy = ((me.y + 1.0) / cell_size) as i32;
        let cc = cell_count as i32;
        let contact = 2.0 * PARTICLE_RADIUS;

        let mut dxa = 0.0;
        let mut dya = 0.0;

        let mut dy = -1i32;
        while dy <= 1 {
            let mut dx = -1i32;
            while dx <= 1 {
                let nx = mcx + dx;
                let ny = mcy + dy;
                if !(nx < 0 || nx >= cc || ny < 0 || ny >= cc) {
                    let cell_idx = (nx as u32 * cell_count + ny as u32) as usize;
                    let n_in_cell = counters[cell_idx];
                    let n = if n_in_cell > GRID_DEPTH as u32 {
                        GRID_DEPTH as u32
                    } else {
                        n_in_cell
                    };
                    let mut slot = 0u32;
                    while slot < n {
                        let j = cells[cell_idx * GRID_DEPTH + slot as usize] as usize;
                        if j != i {
                            // Safety: j was written by grid_build with j < count.
                            let other = unsafe { *positions.add(j) };
                            let xd = me.x - other.x;
                            let yd = me.y - other.y;
                            let d2 = xd * xd + yd * yd;
                            if d2 <= contact * contact {
                                let inv = fast_rsqrt(d2);
                                let d = d2 * inv; // sqrt via x * (1/sqrt(x))
                                let overlap_distance = contact - d;
                                let (nrm_x, nrm_y) = if d > 0.0 {
                                    (xd * inv, yd * inv)
                                } else {
                                    // Index-parity split keeps neighbouring colocated pairs
                                    // from all pushing the same way.
                                    if i & 1 == 0 { (1.0, 0.0) } else { (-1.0, 0.0) }
                                };
                                let correction = RESTITUTION * overlap_distance * 0.25;
                                dxa += correction * nrm_x;
                                dya += correction * nrm_y;
                            }
                        }
                        slot += 1;
                    }
                }
                dx += 1;
            }
            dy += 1;
        }
        // Cap the per-substep shove at MAX_OVERLAP_PUSH; otherwise a deep overlap encodes
        // a huge implicit velocity into (x - ox) and the next verlet call "explodes".
        let mag2 = dxa * dxa + dya * dya;
        let pre_cap_mag = fast_sqrt(mag2);
        diag_max_mag(diag, DIAG_OVL_MAG, pre_cap_mag);
        let cap2 = MAX_OVERLAP_PUSH * MAX_OVERLAP_PUSH;
        if mag2 > cap2 {
            diag_incr(diag, DIAG_OVL_FIRES);
            let scale = MAX_OVERLAP_PUSH / pre_cap_mag;
            dxa *= scale;
            dya *= scale;
        }
        // Safety: i < count; each thread writes only to its own slot, so no aliasing hazard.
        unsafe {
            (*positions.add(i)).x += dxa;
            (*positions.add(i)).y += dya;
        }
    }
}
