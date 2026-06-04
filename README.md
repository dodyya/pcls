# pcls

![Demo gif](https://github.com/dodyya/pcls/releases/download/demo-assets/demo.gif)

A GPU particle simulation. Tens of thousands of particles run as CUDA kernels and render through OpenGL via CUDA-GL interop. Each particle has a hue. Similar hues attract, opposite hues repel.

## Run

You need an NVIDIA GPU, a CUDA toolkit install, and [cuda-oxide](https://github.com/NVlabs/cuda-oxide). Then:

```bash
cargo oxide run --release
```

The target GPU arch is set to `sm_89` in `.cargo/config.toml`. Change it if your card is different.

## Controls

- Left click (or drag): spawn particles
- Right click: clear
- Space: toggle gravity (radial vs. downward)
- M: toggle the hue force
- D: toggle the donut boundary
- S: halt all motion
- Q: quit

## Internals

CUDA kernels handle force evaluation, the uniform spatial grid, and Verlet integration. The grid is updated incrementally rather than rebuilt every step. Particle positions live in a VBO that CUDA writes to directly, so frames never round-trip through CPU memory.
