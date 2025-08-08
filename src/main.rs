use winit::{
    event::{ElementState, Event, MouseButton, WindowEvent},
    event_loop::ControlFlow,
};

mod gfx;
mod grid;
mod gridphx;
mod particles;
mod sweep;
mod sweepphx;
use gfx::Gfx;
use rand::Rng;
use sweepphx::Phx;
const PARTICLE_SIZE: f32 = 0.01;
const WINDOW_SIZE: u32 = 700;

fn main() {
    let (mut gfx, event_loop) = Gfx::new(WINDOW_SIZE, WINDOW_SIZE);

    let mut simulation = Phx::new(PARTICLE_SIZE * 2.0, vec![]);

    let mut last_cursor_pos: Option<(f32, f32)> = None;
    let mut mouse_down = false;
    let mut rng = rand::thread_rng();

    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Poll;

        match &event {
            Event::WindowEvent { window_id, event } => {
                if *window_id != gfx.window.id() {
                    return;
                }
                match event {
                    WindowEvent::CloseRequested => {
                        *control_flow = ControlFlow::Exit;
                        return;
                    }
                    WindowEvent::CursorMoved { position, .. } => {
                        last_cursor_pos = Some((position.x as f32, position.y as f32));
                    }
                    WindowEvent::MouseInput { state, button, .. } => {
                        if *button == MouseButton::Left {
                            mouse_down = *state == ElementState::Pressed;
                        } else if *button == MouseButton::Right && *state == ElementState::Pressed {
                            simulation.clear();
                        }
                    }
                    _ => {}
                }
            }
            Event::RedrawRequested(window_id) => {
                if *window_id != gfx.window.id() {
                    return;
                }

                if mouse_down {
                    if let Some((cursor_x, cursor_y)) = last_cursor_pos {
                        if cursor_x >= 0.0
                            && cursor_x <= WINDOW_SIZE as f32
                            && cursor_y >= 0.0
                            && cursor_y <= WINDOW_SIZE as f32
                        {
                            let sim_x = (cursor_x / WINDOW_SIZE as f32) * 2.0 - 1.0;
                            let sim_y = 1.0 - (cursor_y / WINDOW_SIZE as f32) * 2.0;

                            let vel_x = rng.gen_range(-0.01..0.01);
                            let vel_y = rng.gen_range(-0.01..0.01);

                            simulation.add_particle(sim_x, sim_y, PARTICLE_SIZE, vel_x, vel_y, 1.0);
                        }
                    }
                }

                simulation.step();

                let particles = simulation.get_drawable_particles();
                gfx.clear_frame();
                gfx.draw_particles(&particles);
                gfx.render();
                gfx.request_redraw();
            }
            _ => {}
        }
    });
}
