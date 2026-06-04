use crate::device_view::DeviceView;
use crate::gl_interop::ParticleVbo;
use crate::kernels::kernels as gpu;
use crate::particles::{Particles, RawParticle, XyHue};
use crate::types::Real;
use cuda_core::error::DriverError;
use cuda_core::{CudaContext, CudaStream, DeviceBuffer, LaunchConfig};
use std::sync::Arc;

pub const PARTICLE_RADIUS: Real = 0.0021;

const DT: Real = 1.0 / 60.0;
const SUBSTEPS: usize = 12;
pub const GRAVITY: Real = 1.0;

pub const ANTI_BHOLE: Real = 0.5; // Avoid black holes at center in donut mode
pub const RESTITUTION: Real = 0.9999994;
pub const MAX_V: Real = 0.0015;
pub const MAX_OVERLAP_PUSH: Real = 2.0 * PARTICLE_RADIUS;
pub const GRID_DEPTH: usize = 6;
pub const VELOCITY_DAMPING: Real = 0.999999;
pub const K: Real = 200. * 256. * PARTICLE_RADIUS;
pub const HUE_FORCE_RADIUS: i32 = 7;
pub const HOLE_SIZE: Real = 0.05;
pub const DONUT_SIZE: Real = 1.0;

// Diagnostic slots in the shared 8-slot u32 buffer. Bit-pattern atomic-max on
// non-negative f32 (monotonic) is used for the *_MAG slots; *_FIRES slots are
// straight u32 counters. Single source of truth shared with the kernel writers.
pub const DIAG_OVL_MAG: usize = 0;
pub const DIAG_OVL_FIRES: usize = 1;
pub const DIAG_HUE_MAG: usize = 2;
pub const DIAG_VRP_DISP: usize = 3;
pub const DIAG_VRP_FIRES: usize = 4;
pub const DIAG_VRP_VEL: usize = 5;
pub const DIAG_CON_MAG: usize = 6;
pub const DIAG_CON_FIRES: usize = 7;
pub const DIAG_SLOTS: usize = 8;

pub struct Simulation {
    pub pcls: Particles,
    // Held to keep the CUDA context alive for the lifetime of all device buffers.
    pub _ctx: Arc<CudaContext>,
    pub stream: Arc<CudaStream>,
    // Hot positions (x, y, hue): read by every kernel, written by constrain/overlap/verlet.
    // Also the layout the VBO uses (extract_to_vbo just memcpy_dtod's this into the VBO).
    pub device_pos: DeviceBuffer<XyHue>,
    // Cold integration state (ox, oy, ax, ay): written by gravity/hue_force/verlet,
    // read by verlet only.
    pub device_integ: DeviceBuffer<RawParticle>,
    // Grid: cells holds cell_count^2 * GRID_DEPTH particle ids; counters holds
    // cell_count^2 per-cell slot counts (zeroed before each build).
    pub grid_cells: DeviceBuffer<u32>,
    pub grid_counters: DeviceBuffer<u32>,
    pub cell_count: u32,
    pub module: gpu::LoadedModule,
    // Per-frame diagnostics; slot layout in DIAG_* consts above.
    pub diag: DeviceBuffer<u32>,
    // Lowest host-index whose particle hasn't been pushed to the device yet.
    // add_particle leaves this at the previous count; sync_tail flushes [dirty_start..count].
    dirty_start: usize,
}

impl Simulation {
    pub fn new(cell_size: Real) -> Self {
        let ctx = CudaContext::new(0).expect("CUDA context");
        // cuda-host's embedded loader uses CUDA_OXIDE_TARGET to pick the SM arch
        // when compiling NVVM IR → cubin; sniff the actual device's capability.
        if std::env::var_os("CUDA_OXIDE_TARGET").is_none() {
            let (major, minor) = ctx.compute_capability().expect("compute cap");
            // Safety: process init, single-threaded, no other readers of env yet.
            unsafe {
                std::env::set_var("CUDA_OXIDE_TARGET", format!("sm_{}{}", major, minor));
            }
        }
        let stream = ctx.default_stream();
        let module = gpu::load(&ctx).expect("load embedded kernels");

        let cell_count = (2.0 / cell_size).ceil() as u32;
        let cells_len = (cell_count as usize).pow(2) * GRID_DEPTH;
        let counters_len = (cell_count as usize).pow(2);
        let grid_cells = DeviceBuffer::zeroed(&stream, cells_len).expect("alloc grid cells");
        let grid_counters =
            DeviceBuffer::zeroed(&stream, counters_len).expect("alloc grid counters");
        let device_pos =
            DeviceBuffer::zeroed(&stream, 64).expect("alloc initial device positions");
        let device_integ =
            DeviceBuffer::zeroed(&stream, 64).expect("alloc initial device integ");
        let diag = DeviceBuffer::zeroed(&stream, DIAG_SLOTS).expect("alloc diag");

        Self {
            pcls: Particles::new(60_000),
            _ctx: ctx,
            stream,
            device_pos,
            device_integ,
            grid_cells,
            grid_counters,
            cell_count,
            module,
            diag,
            dirty_start: 0,
        }
    }

    fn reset_diag(&self) {
        let bytes = self.diag.len() * size_of::<u32>();
        // Safety: `self.diag` owns the device allocation; size is its true byte length.
        unsafe {
            cuda_core::memory::memset_d8_async(
                self.diag.cu_deviceptr(),
                0,
                bytes,
                self.stream.cu_stream(),
            )
            .expect("memset diag");
        }
    }

    fn read_diag(&self) -> [u32; DIAG_SLOTS] {
        let mut buf = [0u32; DIAG_SLOTS];
        self.diag
            .copy_to_host(&self.stream, &mut buf)
            .expect("copy diag to host");
        buf
    }

    // Push freshly-added particles [dirty_start..count] from the host mirror to both
    // device buffers. Runs at the top of step(); typically a no-op since add_particle
    // is rare.
    fn sync_tail(&mut self) -> Result<(), DriverError> {
        if self.dirty_start >= self.pcls.count {
            return Ok(());
        }
        let start = self.dirty_start;
        let n = self.pcls.count - start;
        let pos_bytes = size_of::<XyHue>();
        let integ_bytes = size_of::<RawParticle>();
        let pos_dst = self.device_pos.cu_deviceptr() + (start * pos_bytes) as u64;
        let integ_dst = self.device_integ.cu_deviceptr() + (start * integ_bytes) as u64;
        // Safety: device buffers have cap ≥ count (grow ensures this); host slices
        // [start..start+n] are in bounds; copy sizes match.
        unsafe {
            cuda_core::memory::memcpy_htod_async(
                pos_dst,
                self.pcls.positions.as_ptr().add(start),
                n * pos_bytes,
                self.stream.cu_stream(),
            )?;
            cuda_core::memory::memcpy_htod_async(
                integ_dst,
                self.pcls.integs.as_ptr().add(start),
                n * integ_bytes,
                self.stream.cu_stream(),
            )?;
        }
        self.stream.synchronize()?;
        self.dirty_start = self.pcls.count;
        Ok(())
    }

    // Grow the sim's device buffers to at least `new_cap`. Preserves [0..pcls.count]
    // by pulling the device-synced prefix back into the host mirrors, then re-uploading.
    //
    // Only [0..dirty_start] is dtoh'd: anything past dirty_start on the host is freshly
    // added by add_particle and hasn't been uploaded yet, so dtoh would clobber real host
    // data with stale device padding (zero positions, hue 0 = red).
    fn grow_device(&mut self, new_cap: usize) {
        let synced_len = self.dirty_start.min(self.device_pos.len());
        if synced_len > 0 {
            // Safety: device buffers have cap ≥ synced_len; host vectors already cover
            // [0..count] ≥ [0..synced_len] (push() grew them), so the writes are in bounds.
            unsafe {
                cuda_core::memory::memcpy_dtoh_async(
                    self.pcls.positions.as_mut_ptr(),
                    self.device_pos.cu_deviceptr(),
                    synced_len * size_of::<XyHue>(),
                    self.stream.cu_stream(),
                )
                .expect("memcpy positions dtoh on grow");
                cuda_core::memory::memcpy_dtoh_async(
                    self.pcls.integs.as_mut_ptr(),
                    self.device_integ.cu_deviceptr(),
                    synced_len * size_of::<RawParticle>(),
                    self.stream.cu_stream(),
                )
                .expect("memcpy integs dtoh on grow");
            }
            self.stream.synchronize().expect("sync after grow dtoh");
        }
        let mut padded_pos = self.pcls.positions.clone();
        padded_pos.resize(new_cap, XyHue::default());
        let mut padded_integ = self.pcls.integs.clone();
        padded_integ.resize(new_cap, RawParticle::default());
        self.device_pos =
            DeviceBuffer::from_host(&self.stream, &padded_pos).expect("grow device_pos");
        self.device_integ =
            DeviceBuffer::from_host(&self.stream, &padded_integ).expect("grow device_integ");
        self.dirty_start = self.pcls.count;
    }

    fn reset_counters(&self) {
        let bytes = self.grid_counters.len() * size_of::<u32>();
        // Safety: `self.grid_counters` owns the device allocation; size matches.
        unsafe {
            cuda_core::memory::memset_d8_async(
                self.grid_counters.cu_deviceptr(),
                0,
                bytes,
                self.stream.cu_stream(),
            )
            .expect("memset grid counters");
        }
    }

    fn pos_view(&self, n: usize) -> DeviceView<'_, XyHue> {
        DeviceView::from_buffer(&self.device_pos, n)
    }

    fn integ_view(&self, n: usize) -> DeviceView<'_, RawParticle> {
        DeviceView::from_buffer(&self.device_integ, n)
    }

    fn build_grid(&mut self) {
        let n = self.pcls.count;
        if n == 0 {
            return;
        }
        self.reset_counters();
        let pos_view = self.pos_view(n);
        // Safety: cells is written at slot indices given by the atomic counter — disjoint per
        // thread; sizing matches cell_count^2 * GRID_DEPTH.
        unsafe {
            self.module
                .grid_build(
                    &self.stream,
                    LaunchConfig::for_num_elems(n as u32),
                    &pos_view,
                    self.grid_cells.cu_deviceptr() as *mut u32,
                    &self.grid_counters,
                    self.cell_count,
                )
                .expect("grid_build kernel");
        }
    }

    fn launch_gravity(&mut self) {
        let n = self.pcls.count;
        if n == 0 {
            return;
        }
        let g = self.pcls.g_toward_center as u32;
        let pos_view = self.pos_view(n);
        let mut integ_view = self.integ_view(n);
        self.module
            .gravity(
                &self.stream,
                LaunchConfig::for_num_elems(n as u32),
                &pos_view,
                &mut integ_view,
                g,
            )
            .expect("gravity kernel");
    }

    fn launch_constrain(&mut self) {
        let n = self.pcls.count;
        if n == 0 {
            return;
        }
        let d = self.pcls.donut_enabled as u32;
        let mut pos_view = self.pos_view(n);
        self.module
            .constrain(
                &self.stream,
                LaunchConfig::for_num_elems(n as u32),
                &mut pos_view,
                d,
                &self.diag,
            )
            .expect("constrain kernel");
    }

    fn launch_verlet(&mut self, dt: Real) {
        let n = self.pcls.count;
        if n == 0 {
            return;
        }
        let mut pos_view = self.pos_view(n);
        let mut integ_view = self.integ_view(n);
        self.module
            .verlet(
                &self.stream,
                LaunchConfig::for_num_elems(n as u32),
                &mut pos_view,
                &mut integ_view,
                dt,
                &self.diag,
            )
            .expect("verlet kernel");
    }

    fn launch_hue_force(&mut self) {
        let n = self.pcls.count;
        if n == 0 {
            return;
        }
        let pos_view = self.pos_view(n);
        let mut integ_view = self.integ_view(n);
        self.module
            .hue_force(
                &self.stream,
                LaunchConfig::for_num_elems(n as u32),
                &pos_view,
                &mut integ_view,
                &self.grid_cells,
                &self.grid_counters,
                self.cell_count,
                &self.diag,
            )
            .expect("hue_force kernel");
    }

    fn launch_overlap(&mut self) {
        let n = self.pcls.count as u32;
        if n == 0 {
            return;
        }
        // Safety: positions aliases reads of neighbours' x/y with writes to own x/y;
        // the race is intentional iterative refinement — see kernel comment.
        unsafe {
            self.module
                .overlap(
                    &self.stream,
                    LaunchConfig::for_num_elems(n),
                    self.device_pos.cu_deviceptr() as *mut XyHue,
                    &self.grid_cells,
                    &self.grid_counters,
                    self.cell_count,
                    n,
                    &self.diag,
                )
                .expect("overlap kernel");
        }
    }

    pub fn step(&mut self) {
        // PCLS_DEBUG path syncs after each launch so per-kernel timings are meaningful;
        // prod path enqueues the whole substep async.
        use std::sync::atomic::{AtomicU32, AtomicU64, Ordering::Relaxed as R};
        use std::time::Instant;
        static FRAMES: AtomicU64 = AtomicU64::new(0);
        static T_GRAV: AtomicU64 = AtomicU64::new(0);
        static T_HUE: AtomicU64 = AtomicU64::new(0);
        static T_CON: AtomicU64 = AtomicU64::new(0);
        static T_GRID: AtomicU64 = AtomicU64::new(0);
        static T_OVL: AtomicU64 = AtomicU64::new(0);
        static T_VER: AtomicU64 = AtomicU64::new(0);
        static D_OVL_MAG: AtomicU32 = AtomicU32::new(0);
        static D_OVL_FIRES: AtomicU64 = AtomicU64::new(0);
        static D_HUE_MAG: AtomicU32 = AtomicU32::new(0);
        static D_VRP_DISP: AtomicU32 = AtomicU32::new(0);
        static D_VRP_FIRES: AtomicU64 = AtomicU64::new(0);
        static D_VRP_VEL: AtomicU32 = AtomicU32::new(0);
        static D_CON_MAG: AtomicU32 = AtomicU32::new(0);
        static D_CON_FIRES: AtomicU64 = AtomicU64::new(0);

        static DEBUG: std::sync::LazyLock<bool> =
            std::sync::LazyLock::new(|| std::env::var_os("PCLS_DEBUG").is_some());

        self.sync_tail().expect("htod tail");
        self.reset_diag();

        let hue_on = self.pcls.hue_force_enabled;
        let debug = *DEBUG;

        for i in 0..SUBSTEPS {
            if debug {
                let sync = |s: &CudaStream| s.synchronize().expect("phase sync");
                let t_pre_grid = Instant::now();
                self.build_grid();
                sync(&self.stream);
                let t_post_grid = Instant::now();

                self.launch_gravity();
                sync(&self.stream);
                let t1 = Instant::now();
                if hue_on && i % 4 == 0 {
                    self.launch_hue_force();
                    sync(&self.stream);
                }
                let t2 = Instant::now();
                self.launch_constrain();
                sync(&self.stream);
                let t3 = Instant::now();

                self.build_grid();
                sync(&self.stream);
                let t4 = Instant::now();

                self.launch_overlap();
                sync(&self.stream);
                let t5 = Instant::now();
                self.launch_verlet(DT / SUBSTEPS as Real);
                sync(&self.stream);
                let t6 = Instant::now();

                T_GRAV.fetch_add((t1 - t_post_grid).as_nanos() as u64, R);
                T_HUE.fetch_add((t2 - t1).as_nanos() as u64, R);
                T_CON.fetch_add((t3 - t2).as_nanos() as u64, R);
                T_GRID.fetch_add(
                    ((t_post_grid - t_pre_grid) + (t4 - t3)).as_nanos() as u64,
                    R,
                );
                T_OVL.fetch_add((t5 - t4).as_nanos() as u64, R);
                T_VER.fetch_add((t6 - t5).as_nanos() as u64, R);
            } else {
                self.build_grid();
                self.launch_gravity();
                if hue_on && i % 4 == 0 {
                    self.launch_hue_force();
                }
                self.launch_constrain();
                self.build_grid();
                self.launch_overlap();
                self.launch_verlet(DT / SUBSTEPS as Real);
            }
        }

        let d = self.read_diag();
        // Bit-pattern atomic-max is monotonic for non-negative f32, which all magnitude slots are.
        D_OVL_MAG.fetch_max(d[DIAG_OVL_MAG], R);
        D_OVL_FIRES.fetch_add(d[DIAG_OVL_FIRES] as u64, R);
        D_HUE_MAG.fetch_max(d[DIAG_HUE_MAG], R);
        D_VRP_DISP.fetch_max(d[DIAG_VRP_DISP], R);
        D_VRP_FIRES.fetch_add(d[DIAG_VRP_FIRES] as u64, R);
        D_VRP_VEL.fetch_max(d[DIAG_VRP_VEL], R);
        D_CON_MAG.fetch_max(d[DIAG_CON_MAG], R);
        D_CON_FIRES.fetch_add(d[DIAG_CON_FIRES] as u64, R);

        let n = FRAMES.fetch_add(1, R) + 1;
        if debug && n.is_multiple_of(60) {
            static LAST_PRINT: std::sync::Mutex<Option<Instant>> = std::sync::Mutex::new(None);
            static P_GRAV: AtomicU64 = AtomicU64::new(0);
            static P_HUE: AtomicU64 = AtomicU64::new(0);
            static P_CON: AtomicU64 = AtomicU64::new(0);
            static P_GRID: AtomicU64 = AtomicU64::new(0);
            static P_OVL: AtomicU64 = AtomicU64::new(0);
            static P_VER: AtomicU64 = AtomicU64::new(0);

            let now = Instant::now();
            let fps = {
                let mut lp = LAST_PRINT.lock().unwrap();
                let fps = lp.map_or(0.0, |prev| 60.0 / (now - prev).as_secs_f64());
                *lp = Some(now);
                fps
            };
            let win = (60 * SUBSTEPS as u64) as f64 * 1000.0;
            let delta = |prev: &AtomicU64, cur: &AtomicU64| {
                let now = cur.load(R);
                let was = prev.swap(now, R);
                (now - was) as f64 / win
            };
            eprintln!(
                "[sim n={:>5} pcls f={:>4}] fps={:>5.1} grav={:>5.1} hue={:>6.1} con={:>5.1} grid={:>5.1} ovl={:>6.1} ver={:>5.1}   (µs/substep, last 60 frames, hue on={})",
                self.pcls.count,
                n,
                fps,
                delta(&P_GRAV, &T_GRAV),
                delta(&P_HUE, &T_HUE),
                delta(&P_CON, &T_CON),
                delta(&P_GRID, &T_GRID),
                delta(&P_OVL, &T_OVL),
                delta(&P_VER, &T_VER),
                self.pcls.hue_force_enabled,
            );
            let dec = |a: &AtomicU32| f32::from_bits(a.swap(0, R));
            let ovl_fires = D_OVL_FIRES.swap(0, R);
            let vrp_fires = D_VRP_FIRES.swap(0, R);
            let con_fires = D_CON_FIRES.swap(0, R);
            eprintln!(
                "[diag] overlap_pre_cap_max={:.5} overlap_cap_fires={} hue_force_max={:.5} verlet_pre_clamp_max={:.5} verlet_clamp_fires={} carried_v_max={:.5} constrain_shove_max={:.5} constrain_big_shoves={}    (caps MAX_OVERLAP_PUSH={:.5} MAX_V={:.5})",
                dec(&D_OVL_MAG),
                ovl_fires,
                dec(&D_HUE_MAG),
                dec(&D_VRP_DISP),
                vrp_fires,
                dec(&D_VRP_VEL),
                dec(&D_CON_MAG),
                con_fires,
                MAX_OVERLAP_PUSH,
                MAX_V,
            );
        }
    }

    // Copy the current device positions into the GL-mapped VBO. Positions are already
    // XyHue-shaped (same layout the VAO binds to), so this is a single device-to-device
    // memcpy — no per-thread kernel needed. Caller must grow the VBO to fit pcls.count.
    pub fn extract_to_vbo(&mut self, vbo: &mut ParticleVbo) -> Result<(), DriverError> {
        let n = self.pcls.count;
        if n == 0 {
            return Ok(());
        }
        assert!(n <= vbo.cap);
        let dst = vbo.map(&self.stream)?;
        let bytes = n * size_of::<XyHue>();
        // Safety: `dst` is the VBO's CUDA-mapped pointer (valid until vbo.unmap), the
        // source is our typed positions buffer with len ≥ n, both in the same context.
        unsafe {
            cuda_core::memory::memcpy_dtod_async(
                dst,
                self.device_pos.cu_deviceptr(),
                bytes,
                self.stream.cu_stream(),
            )?;
        }
        // Sync before unmap so GL sees CUDA's writes when it next binds the VBO.
        self.stream.synchronize()?;
        vbo.unmap(&self.stream)?;
        Ok(())
    }

    pub fn is_hue_force_enabled(&self) -> bool {
        self.pcls.hue_force_enabled
    }

    pub fn add_particle(&mut self, x: Real, y: Real, hue: Real) {
        self.pcls.push(x, y, hue);
        if self.pcls.count > self.device_pos.len() {
            self.grow_device(self.pcls.count.next_power_of_two().max(64));
        }
    }

    pub fn clear(&mut self) {
        self.pcls.clear();
        self.dirty_start = 0;
    }

    pub fn toggle_gravity(&mut self) {
        self.pcls.g_toward_center = !self.pcls.g_toward_center;
    }

    pub fn toggle_hue_force(&mut self) {
        self.pcls.hue_force_enabled = !self.pcls.hue_force_enabled;
    }

    pub fn toggle_donut(&mut self) {
        self.pcls.donut_enabled = !self.pcls.donut_enabled;
    }

    pub fn stop(&mut self) {
        let n = self.pcls.count;
        if n == 0 {
            return;
        }
        let pos_view = self.pos_view(n);
        let mut integ_view = self.integ_view(n);
        self.module
            .zero_velocity(
                &self.stream,
                LaunchConfig::for_num_elems(n as u32),
                &pos_view,
                &mut integ_view,
            )
            .expect("zero_velocity kernel");
    }
}
