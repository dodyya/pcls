# PCLS - Particle Physics Simulation

![Demo gif](quickerdemo.gif)
![Demo video](demo.mp4)
A decnetly fast, multithreaded particle physics simulation engine written in Rust, capable of simulating tens of thousands of particles in real-time with gravity, collision detection, and other miscellaneous forces. (GIF displays slower in GitHub than it should)

## Overview

PCLS (Particles) is a simple physics simulation that implements Verlet integration for stable particle motion, spatial partitioning using uniform grids for O(n) collision detection, multithreaded collision resolution with atomic data types for memory safety, Coulomb electrostatic forces between charged particles, variations in mass, radius, and charge, along with a responsive, punchy visualization.

The simulation can handle about 50,000 particles at 60 FPS on my machine with Coulomb turned off.

## Controls

Run using Cargo: `cargo run --release`.
Left click to spawn particles continuously, right click to clear screen of particles. Toggle gravity between radially towards the center of the screen and downwards with SPACE.
Halt all particles with S. Toggle particle constraint with D, turn electromagnetic effects off with M.

## Technical Implementation

Particles are stored as an array of structs, storing current position, previous frame position, and acceleration.
Each subframe, forces (gravity, Coulomb) determine the acceleration of the particle. Collision detection and resolution forces particles apart, and constraint resolution forces particles to remain within bounds. Verlet integration is performed to determine each particle's new position.
Due to the way Verlet integration works, the velocity of a particle is computed dynamically, and does not need to be stored. This means that less invariants have to be upheld -- it is easier to make sure that positions are valid and collisions have happened if you don't have to make sure velocity is up-to-date and valid. I'm saying this from experience -- I started this project in February of last year with no resources or research into physics simulations, and not knowing about Verlet meant I struggled with invalid behaviors. I resorted to heavily clamping velocities, splitting collisions into positional and velocity resolution, and dealing with the complex interactions of particle positions and velocities on contact. Constraining particles to any sort of non-rectangle shape also needed to be done manually. Pivoting to Verlet integration (after stumbling upon a YouTube video describing a similar project, in C++) fixed every single one of these issues immediately.
As long as each step did a reasonable job moving the particle from an invalid state towards a valid one, the Verlet step would continue to do so. Overlaps happened in 1 iteration. Velocity did not need to be aggressively clamped. Any target constraint could be achieved, so long as there could be a function to take in an out-of-bounds particle and move it to the closest in-bounds point. Furthermore, particle 'crushing' towards the bottom of the screen while particles were at rest did not take place, as with a previous-position essentially identical to the current position, and with micro-collisions and substeps avoiding excessive overlaps, the simulation became stable.

## Optimizations

The Array3D implementation provides superior cache locality compared to nested vectors. Rather than storing data as Vec<Vec<Vec<T>>> which would scatter memory allocations across the heap, my custom Array3D allocates a single contiguous block and uses index arithmetic to access elements. This ensures that neighboring grid cells are physically adjacent in memory, dramatically improving cache hit rates during neighbor searches. It should also be clear that particle structs themselves were never copied or moved, and their indices in the Particles array were instead used.

The particle system uses Array-of-Structures (AoS) layout where each particle contains all its properties (position, velocity, acceleration, mass, charge) together in memory. While I did initially turn to SoA because "it should just be faster", I never utilized the minor SIMD benefits of it. Benchmarking with cargo-flamegraph, I found that non-SIMD-able steps of the simulation were by far the most expensive, partially because they required many fields of the same particle to be accessed at once. AoS could mean that more of these fields would be more likely to experience cache hits. Testing AoS vs SoA showed a ~20% frame rate increase, and I pivoted to Array of Structs.

## Atomics and Threads

Having had minimal experience with parallelism and concurrency prior to this project, and reaching that stage of the project on a plane with only the Rust book and a locally hosted model to help me, I attempted to derive a simple multithreading solution. Arc<T> allowed me to access data between threads, and RwLock<T> allowed me to safely share data between threads. However, even though I never struggled with conflicting locks on data, the sheer amount of lock operations brought performance down. Switching to Atomics, particularly AtomicUsize and atomic float(s) from a crate of the same name meant that I could avoid the overhead of locks altogether. Furthermore, going from a SeqCst (guaranteeing threads see all sequentially consistent operations in the same order) to a Relaxed (not) ordering of atomic operations sped up performance further.

## Graphics Implementation

The simulation uses the pixels crate for simple, efficient graphics rendering. Rather than using complex graphics APIs like OpenGL or Vulkan, pixels provides direct access to a pixel buffer that maps to the display surface. This approach offers several advantages for this application.

The gfx.rs module implements basic drawing primitives (copied from a previous graphics project) directly by manipulating pixel values in the buffer. Circle rendering uses Bresenham's circle algorithm, and colors are found in a lookup table.

The pixels crate handles the platform-specific details of presenting the pixel buffer to the screen while providing a simple, consistent interface with winit. This allowed me to focus on the physics calculations rather than graphics API complexity. The direct pixel manipulation approach also enables easy frame capture for animation recording by simply copying the pixel buffer data, which I used to test another image pipeline project (qoi). Unstable features related to graphics: press P to dump a single frame's dimension-prefixed display buffer representation of the visible window to stdout, and press R to toggle "recording" to stdout. Example: `./pcls | ./qoi write -fn out.qoi` pipes the pixel buffer data to the qoi encoder, which writes each frame to an individual numbered .qoi file, which can then be processed into a video/gif.

## Dependencies

The simulation relies on `pixels` for hardware-accelerated pixel buffer rendering, `winit` for cross-platform windowing and event handling, `rayon` for selective data parallelism in specific hot paths, `atomic_float` for lock-free atomic floating point operations essential for thread safety, and `rand` for pseudo-random number generation during particle placement.

## Usage

Run the simulation with `cargo run --release` for passable performance. The release mode is essential as debug builds cannot achieve real-time framerates with thousands of particles.

For animation creation, enable recording mode by pressing R during simulation, then convert the exported frame sequence to video using ffmpeg with commands like `ffmpeg -framerate 60 -i frames%05d.png -pix_fmt yuv420p output.mp4`.

## Future Enhancements

The Barnes-Hut algorithm could provide O(n log n) complexity for long-range forces using hierarchical spatial subdivision. A QuadTree structure with center of mass calculations, total mass tracking, spatial bounds, optional children nodes, and particle lists would enable adaptive precision based on distance and massive scalability improvements for 100k+ particles.

GPU acceleration through WGPU compute shaders could enable parallel force calculation on thousands of GPU cores, shared memory optimization for collision detection, and async CPU-GPU pipelining for maximum throughput. Advanced physics implementations could include fluid dynamics with pressure and viscosity, soft body simulation with spring constraints, heat transfer and thermal dynamics, and magnetic field interactions with Lorentz forces.

Using more unsafe code and accessing data directly (raw pointers instead of Arc) could increase performance further, but I had admittedly reached a point of "good enough" where it could've made more sense to optimize my graphics rendering pipeline first, and I really did not care to move the focus to that aspect of the project.

## Sources/inspiration

Verlet suggestion:
https://www.youtube.com/watch?v=lS_qeBy3aQI

Multithreading:
https://www.youtube.com/watch?v=9IULfQH7E90

P.S. I promise I did not simply rewrite the code of the individual in the video.
