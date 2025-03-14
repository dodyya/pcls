#![allow(unused)]
use winit::{
    event::{ElementState, Event, MouseButton, WindowEvent},
    event_loop::ControlFlow,
};

mod gfx;
mod gridphx;
mod particles;
mod sweepphx;
use gfx::Gfx;
use gridphx::Phx as GridPhx;
use rand::Rng;
use sweepphx::Phx as SweepPhx;

const PARTICLE_SIZE: f32 = 0.01;
const WINDOW_SIZE: u32 = 700; // Smaller window size to fit two windows
const WINDOW_SPACING: i32 = 20; // Spacing between windows

fn main() {
    // Create event loop
    let event_loop = winit::event_loop::EventLoop::new();

    // Create two windows
    let window1 = winit::window::WindowBuilder::new()
        .with_title("Sweep Physics")
        .with_inner_size(winit::dpi::PhysicalSize::new(WINDOW_SIZE, WINDOW_SIZE))
        .with_position(winit::dpi::PhysicalPosition::new(50, 100))
        .with_resizable(false)
        .build(&event_loop)
        .unwrap();

    let window2 = winit::window::WindowBuilder::new()
        .with_title("Grid Physics")
        .with_inner_size(winit::dpi::PhysicalSize::new(WINDOW_SIZE, WINDOW_SIZE))
        .with_position(winit::dpi::PhysicalPosition::new(
            50 + WINDOW_SIZE as i32 + WINDOW_SPACING,
            100,
        ))
        .with_resizable(false)
        .build(&event_loop)
        .unwrap();

    // Create pixel buffers for the windows
    let pixels1 = pixels::Pixels::new(
        WINDOW_SIZE,
        WINDOW_SIZE,
        pixels::SurfaceTexture::new(WINDOW_SIZE, WINDOW_SIZE, &window1),
    )
    .unwrap();

    let pixels2 = pixels::Pixels::new(
        WINDOW_SIZE,
        WINDOW_SIZE,
        pixels::SurfaceTexture::new(WINDOW_SIZE, WINDOW_SIZE, &window2),
    )
    .unwrap();

    // We need to create separate event loops for each Gfx instance
    // But winit only allows one event loop per application
    // So we'll modify our approach to use a modified version of Gfx

    // Create a modified Gfx struct just for this example
    struct GfxWithWindow {
        window: winit::window::Window,
        pixels: pixels::Pixels,
        width: u32,
        height: u32,
    }

    impl GfxWithWindow {
        fn draw_particles(&mut self, particles: &[(f32, f32, f32)]) {
            for (x, y, r) in particles {
                gfx::fill_circle(self.pixels.frame_mut(), self.width, (*x, *y), *r);
            }
        }

        fn clear_frame(&mut self) {
            let black = [0, 0, 0, 255].repeat(self.pixels.frame().len() / 4);
            self.pixels.frame_mut().copy_from_slice(&black);
        }

        fn render(&mut self) {
            self.pixels.render().unwrap();
        }

        fn request_redraw(&self) {
            self.window.request_redraw();
        }
    }

    // Create our custom Gfx instances
    let mut gfx1 = GfxWithWindow {
        window: window1,
        pixels: pixels1,
        width: WINDOW_SIZE,
        height: WINDOW_SIZE,
    };

    let mut gfx2 = GfxWithWindow {
        window: window2,
        pixels: pixels2,
        width: WINDOW_SIZE,
        height: WINDOW_SIZE,
    };

    // Create two simulation instances
    let mut sweep_simulation = SweepPhx::new(PARTICLE_SIZE * 2.0, vec![]);
    let mut grid_simulation = GridPhx::new(PARTICLE_SIZE * 4.0, vec![]);

    // We'll track the last known cursor position and whether the mouse is pressed.
    let mut last_cursor_pos: Option<(f32, f32)> = None;
    let mut mouse_down = false;
    let mut active_window: Option<u32> = None; // Track which window is receiving input
    let mut rng = rand::thread_rng();

    event_loop.run(move |event, target, control_flow| {
        *control_flow = ControlFlow::Poll;

        match &event {
            Event::WindowEvent { window_id, event } => {
                let is_window1 = window_id == &gfx1.window.id();
                let is_window2 = window_id == &gfx2.window.id();

                match event {
                    WindowEvent::CloseRequested => {
                        *control_flow = ControlFlow::Exit;
                        return;
                    }
                    // Update the last cursor position and track which window is active
                    WindowEvent::CursorMoved { position, .. } => {
                        if is_window1 || is_window2 {
                            last_cursor_pos = Some((position.x as f32, position.y as f32));
                            active_window = Some(if is_window1 { 1 } else { 2 });
                        }
                    }
                    // When the left or right mouse button is pressed or released, update our flag.
                    WindowEvent::MouseInput { state, button, .. } => {
                        if *button == MouseButton::Left {
                            mouse_down = *state == ElementState::Pressed;
                            if mouse_down {
                                active_window = Some(if is_window1 { 1 } else { 2 });
                            }
                        } else if *button == MouseButton::Right && *state == ElementState::Pressed {
                            // Clear all particles in BOTH windows when right-clicking either window
                            sweep_simulation.particles.clear();
                            sweep_simulation.broadphase.clear();
                            grid_simulation.particles.clear();
                            grid_simulation.grid.clear();
                        }
                    }
                    _ => {}
                }
            }
            Event::RedrawRequested(window_id) => {
                let is_window1 = window_id == &gfx1.window.id();
                let is_window2 = window_id == &gfx2.window.id();

                // Add particles if mouse is down and we know which window is active
                if mouse_down && active_window.is_some() {
                    if let Some((cursor_x, cursor_y)) = last_cursor_pos {
                        // Convert window coordinates (0 to window_size) to simulation coordinates (-1 to 1).
                        if cursor_x >= 0.0
                            && cursor_x <= WINDOW_SIZE as f32
                            && cursor_y >= 0.0
                            && cursor_y <= WINDOW_SIZE as f32
                        {
                            let sim_x = (cursor_x / WINDOW_SIZE as f32) * 2.0 - 1.0;
                            let sim_y = 1.0 - (cursor_y / WINDOW_SIZE as f32) * 2.0;

                            // Create particle data tuple with consistent velocities for all simulations
                            let vel_x = rng.gen_range(-0.01..0.01);
                            let vel_y = rng.gen_range(-0.01..0.01);
                            let particle_data = (
                                sim_x,         // center_x
                                sim_y,         // center_y
                                PARTICLE_SIZE, // radius
                                vel_x,         // velocity_x
                                vel_y,         // velocity_y
                                1.0,           // mass
                            );

                            // Add particle to both simulations at the same position
                            // First add to sweep simulation
                            let new_index_sweep = sweep_simulation.particles.count;
                            sweep_simulation.particles.center_x.push(particle_data.0);
                            sweep_simulation.particles.center_y.push(particle_data.1);
                            sweep_simulation.particles.radius.push(particle_data.2);
                            sweep_simulation.particles.velocity_x.push(particle_data.3);
                            sweep_simulation.particles.velocity_y.push(particle_data.4);
                            sweep_simulation.particles.mass.push(particle_data.5);
                            sweep_simulation.particles.count += 1;
                            sweep_simulation
                                .broadphase
                                .rebuild(&sweep_simulation.particles);

                            // Then add to grid simulation
                            let new_index_grid = grid_simulation.particles.count;
                            grid_simulation.particles.center_x.push(particle_data.0);
                            grid_simulation.particles.center_y.push(particle_data.1);
                            grid_simulation.particles.radius.push(particle_data.2);
                            grid_simulation.particles.velocity_x.push(particle_data.3);
                            grid_simulation.particles.velocity_y.push(particle_data.4);
                            grid_simulation.particles.mass.push(particle_data.5);
                            grid_simulation.particles.count += 1;
                            grid_simulation.grid.insert(
                                new_index_grid,
                                particle_data.0,
                                particle_data.1,
                            );
                        }
                    }
                }

                // Update simulations
                let start_sweep = std::time::Instant::now();
                sweep_simulation.update();
                println!(
                    "Sweep update took: {:?} for {} particles",
                    start_sweep.elapsed(),
                    sweep_simulation.particles.count
                );

                let start_grid = std::time::Instant::now();
                grid_simulation.update();
                println!(
                    "Grid update took: {:?} for {} particles",
                    start_grid.elapsed(),
                    grid_simulation.particles.count
                );

                // Draw particles for each window
                if is_window1 {
                    let particles = sweep_simulation.get_drawable_particles();
                    gfx1.clear_frame();
                    gfx1.draw_particles(&particles);
                    gfx1.render();
                    gfx1.request_redraw();
                }

                if is_window2 {
                    let particles = grid_simulation.get_drawable_particles();
                    gfx2.clear_frame();
                    gfx2.draw_particles(&particles);
                    gfx2.render();
                    gfx2.request_redraw();
                }
            }
            _ => {}
        }
    });
}
