use winit::{
    dpi::PhysicalSize,
    event::{ElementState, Event, MouseButton, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::WindowBuilder,
};

mod gfx;
mod phx;
use gfx::Gfx;
use phx::{Particle, Phx, Vec2};

const PARTICLE_SIZE: f32 = 0.1;

fn main() {
    // Create graphics and simulation.
    let (mut gfx, event_loop) = Gfx::new(1200, 1200);

    // Create some initial particles.
    // Make 100 particles in a loop.
    let mut particles = Vec::new();

    let mut simulation = Phx::new(PARTICLE_SIZE * 2.0, particles);

    // We'll track the last known cursor position.
    let mut last_cursor_pos: Option<(f32, f32)> = None;
    // Our window size in pixels.
    let window_size = 1200.0;

    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Poll;

        match &event {
            Event::WindowEvent { event, .. } => match event {
                WindowEvent::CloseRequested => {
                    *control_flow = ControlFlow::Exit;
                    return;
                }
                // Update last_cursor_pos when the mouse moves.
                WindowEvent::CursorMoved { position, .. } => {
                    last_cursor_pos = Some((position.x as f32, position.y as f32));
                }
                // On left mouse click, add a new particle.
                WindowEvent::MouseInput { state, button, .. } => {
                    if *state == ElementState::Pressed && *button == MouseButton::Left {
                        if let Some((cursor_x, cursor_y)) = last_cursor_pos {
                            // Convert from window coordinates (0 to 1200) to simulation coordinates (-1 to 1).
                            // Note: we flip y because window y=0 is at the top.
                            let sim_x = (cursor_x / window_size) * 2.0 - 1.0;
                            let sim_y = 1.0 - (cursor_y / window_size) * 2.0;

                            let new_particle = Particle {
                                center: Vec2 { x: sim_x, y: sim_y },
                                radius: PARTICLE_SIZE,
                                // For now, new particles start at rest.
                                velocity: Vec2 { x: 0.0, y: 0.0 },
                                mass: 1.0,
                            };

                            // Push the new particle into the simulation.
                            simulation.particles.push(new_particle.clone());
                            // Its index is the last index in the vector.
                            let new_index = simulation.particles.len() - 1;
                            // Insert into the hash grid.
                            simulation.grid.insert(new_index, &new_particle);
                        }
                    }
                }
                _ => {}
            },
            Event::RedrawRequested(_) => {
                // Run simulation update.
                simulation.update();

                // Prepare drawing.
                let particles_draw: Vec<(f32, f32, f32)> = simulation.get_drawable_particles();

                gfx.clear_frame();
                gfx.draw_particles(&particles_draw);
                gfx.render();
            }
            _ => {}
        }

        gfx.request_redraw();
    });
}
