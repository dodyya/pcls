use winit::{
    event::{ElementState, Event, MouseButton, WindowEvent},
    event_loop::ControlFlow,
};
//random

mod gfx;
mod phx;
use gfx::Gfx;
use phx::{Particle, Phx, Vec2};
use rand::Rng;

const PARTICLE_SIZE: f32 = 0.01;

fn main() {
    // Create graphics and simulation.
    let (mut gfx, event_loop) = Gfx::new(1200, 1200);

    let mut simulation = Phx::new(PARTICLE_SIZE * 2.0, vec![]);

    // We'll track the last known cursor position and whether the mouse is pressed.
    let mut last_cursor_pos: Option<(f32, f32)> = None;
    let mut mouse_down = false;
    let window_size = 1200.0;

    let mut rng = rand::thread_rng();

    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Poll;

        match &event {
            Event::WindowEvent { event, .. } => match event {
                WindowEvent::CloseRequested => {
                    *control_flow = ControlFlow::Exit;
                    return;
                }
                // Update the last cursor position.
                WindowEvent::CursorMoved { position, .. } => {
                    last_cursor_pos = Some((position.x as f32, position.y as f32));
                }
                // When the left or right mouse button is pressed or released, update our flag.
                WindowEvent::MouseInput { state, button, .. } => {
                    if *button == MouseButton::Left {
                        mouse_down = *state == ElementState::Pressed;
                    } else if *button == MouseButton::Right && *state == ElementState::Pressed {
                        // Clear all particles on right-click
                        simulation.particles.clear();
                        simulation.grid.clear();
                    }
                }
                _ => {}
            },
            Event::RedrawRequested(_) => {
                // Optionally, add a new particle every frame while dragging.
                if mouse_down {
                    if let Some((cursor_x, cursor_y)) = last_cursor_pos {
                        // Convert window coordinates (0 to 1200) to simulation coordinates (-1 to 1).
                        // If mouse is outside window, return
                        if cursor_x < 0.0
                            || cursor_x > window_size
                            || cursor_y < 0.0
                            || cursor_y > window_size
                        {
                            return;
                        }
                        let sim_x = (cursor_x / window_size) * 2.0 - 1.0;
                        let sim_y = 1.0 - (cursor_y / window_size) * 2.0;

                        let new_particle = Particle {
                            center: Vec2 { x: sim_x, y: sim_y },
                            radius: PARTICLE_SIZE,
                            // Random velocity
                            velocity: Vec2 {
                                x: rng.gen_range(-0.01..0.01),
                                y: rng.gen_range(-0.01..0.01),
                            },
                            mass: 1.0,
                        };

                        simulation.particles.push(new_particle.clone());
                        let new_index = simulation.particles.len() - 1;
                        simulation.grid.insert(new_index, &new_particle);
                    }
                }
                //Begin timing
                // let start = std::time::Instant::now();

                // Run simulation update.
                simulation.update();
                // Print how long update took
                // println!("Update took: {:?}", start.elapsed());

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
