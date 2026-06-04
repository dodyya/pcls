# pcls

![Demo gif](https://github.com/dodyya/pcls/releases/download/demo-assets/demo.gif)

A particle simulation. Tens of thousands of particles run on the CPU across every core and render as OpenGL points. Each particle has a hue. Similar hues attract, opposite hues repel.

There's also a [`cuda`](https://github.com/dodyya/pcls/tree/cuda) branch that runs the same simulation as CUDA kernels on the GPU.

## Run

```bash
cargo run --release
```

## Controls

- Left click (or drag): spawn particles
- Right click: clear
- Space: toggle gravity (radial vs. downward)
- M: toggle the hue force
- D: toggle the donut boundary
- S: halt all motion
- Q: quit

## Internals

Particles update with Verlet integration over a uniform spatial grid for O(n) collision detection. The heavy work is parallelized with `rayon`: overlaps resolve in two passes over alternating grid columns, so neighboring columns never touch the same particle at once, and the hue force accumulates through lock-free atomics. Positions and hues upload to a VBO and draw as GL point sprites, with HSV-to-RGB and the circular mask handled in the fragment shader.
