# PCLS - Particle Physics Simulation

![Demo gif](demo.gif)

PCLS is a fast, multithreaded particle physics simulation engine built in Rust. It simulates tens of thousands of particles in real-time, with gravity, collision detection, electrostatic (Coulomb), and additional forces. The engine uses efficient spatial partitioning and multithreaded collision resolution to enable massive scale, and provides direct, punchy pixel-based rendering for a responsive visual experience.

## Features

- Real-time simulation of 50,000+ particles (depending on hardware)
- Gravity (radial or vertical), Coulomb (electrostatic) forces
- Per-particle mass, radius, charge
- O(n) collision detection using uniform grid spatial partitioning
- Multithreaded force & collision resolution (lock-free where possible)
- Verlet integration for stable, efficient physics
- Responsive pixel graphics via the `pixels` and `winit` crates
- Fast animation/video export support

## Controls

- **Run:** `cargo run --release`
- **Left Click:** Spawn particles continuously
- **Right Click:** Clear all particles
- **Space:** Toggle gravity mode
- **S:** Halt all particles
- **Q:** Quit application
- **D:** Toggle particle constraint
- **M:** Toggle electromagnetic effects
- **R:** Toggle frame recording (export pixel buffer to stdout for use with `qoip`)
- **P:** Dump current frame to stdout

## Technical Overview

- Particles are stored as structs with their state: current/previous positions, acceleration
- Force calculations and constraint resolution each tick
- Collision/constraint enforcement and Verlet integration drive particle state
- Data-oriented design (Array-of-Structures) for cache efficiency
- Uses threads and atomics for parallel performance
- Simple, direct pixel drawing for fast, OS-independent display

## Dependencies

- `pixels` - hardware-accelerated rendering
- `winit` - cross-platform window/events
- `rayon`, `atomic_float` - efficient parallelism, lock-free shared state
- `rand` - for random initialization

## Usage

Run with release optimizations:

```bash
cargo run --release
```

- Press `R` while running to enable frame recording (to stdout)
- Pipe output to `qoip` for video encoding

## Further Development Ideas

- Barnes-Hut/QuadTree for scalable long-range force calculations
- GPU compute (WGPU) for massive parallel acceleration
- More complex physics modes: fluids, soft bodies, thermal/magnetic effects

## References

- Verlet Integration: https://www.youtube.com/watch?v=lS_qeBy3aQI
- Multithreading: https://www.youtube.com/watch?v=9IULfQH7E90
