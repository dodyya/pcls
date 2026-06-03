mod app;
mod array;
mod grid;
mod maybe_id;
mod particles;
mod render;
mod sim;
mod types;

use crate::app::App;
fn main() {
    if std::env::var("BENCH").is_ok() {
        bench();
        return;
    }
    App::run();
}

fn bench() {
    use crate::sim::{Simulation, PARTICLE_RADIUS};
    use crate::types::{Real, TAU};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::time::Instant;

    let n: usize = std::env::var("BENCH_PARTICLES")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(5000);
    let steps: usize = std::env::var("BENCH_STEPS")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(300);
    let seed: u64 = std::env::var("BENCH_SEED")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(42);
    let warmup: usize = std::env::var("BENCH_WARMUP")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(60);

    let mut rng = StdRng::seed_from_u64(seed);
    let mut sim = Simulation::new(PARTICLE_RADIUS * 2.0);

    if std::env::var("BENCH_HUE").is_ok() {
        sim.toggle_hue_force();
    }

    // Default donut mode is on; seed an annulus that fits between hole and rim.
    for _ in 0..n {
        let theta: Real = rng.gen_range(0.0..TAU);
        let r: Real = rng.gen_range(0.32..0.95);
        let x = r * theta.cos();
        let y = r * theta.sin();
        let hue: Real = rng.gen_range(0.0..1.0);
        sim.add_particle(x, y, hue);
    }

    eprintln!("BENCH seeded={} warmup={} steps={} seed={}", n, warmup, steps, seed);

    for _ in 0..warmup { sim.step(); }

    let t0 = Instant::now();
    for _ in 0..steps { sim.step(); }
    let total_ms = t0.elapsed().as_secs_f64() * 1000.0;

    println!(
        "BENCH_TOTAL particles={} steps={} total_ms={:.3} ms_per_step={:.4}",
        n, steps, total_ms, total_ms / steps as f64
    );
}
