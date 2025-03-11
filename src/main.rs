use winit::{
    event::{ElementState, Event, MouseButton, WindowEvent},
    event_loop::ControlFlow,
};
mod gfx;
mod phx;
use gfx::Gfx;
use phx::Phx;
use rand::Rng;
const PARTICLE_SIZE: f32 = 0.01;
fn main() {
    let (mut gfx, event_loop) = Gfx::new(1200, 1200);
    let mut simulation = Phx::new(PARTICLE_SIZE * 2.0);
    let mut last_cursor_pos: Option<(f32, f32)> = None;
    let mut mouse_down = false;
    let window_size = 1200.0;
    let mut rng = rand::thread_rng();
    simulation.push(0.5, 0.0, PARTICLE_SIZE, -0.01, 0.05, 1.0);
    simulation.push(-0.5, 0.0, PARTICLE_SIZE, 0.01, 0.05, 1.0);
    event_loop.run(move |event, _, control_flow| {
        *control_flow = ControlFlow::Poll;
        match &event {
            Event::WindowEvent { event, .. } => match event {
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
                        simulation.particles.clear();
                        simulation.grid.clear();
                    }
                }
                _ => {}
            },
            Event::RedrawRequested(_) => {
                if mouse_down {
                    if let Some((cursor_x, cursor_y)) = last_cursor_pos {
                        if cursor_x < 0.0
                            || cursor_x > window_size
                            || cursor_y < 0.0
                            || cursor_y > window_size
                        {
                            return;
                        }
                        let sim_x = (cursor_x / window_size) * 2.0 - 1.0;
                        let sim_y = 1.0 - (cursor_y / window_size) * 2.0;
                        simulation.push(
                            sim_x,
                            sim_y,
                            PARTICLE_SIZE,
                            rng.gen_range(-0.01..0.01),
                            rng.gen_range(-0.01..0.01),
                            1.0,
                        );
                    }
                }
                simulation.update();
                let particles_draw: Vec<(f32, f32, f32)> = simulation.get_drawable_particles();
                gfx.clear_frame();
                gfx.draw_particles(&particles_draw);
                gfx.render();
                //sleep for 500ms
                // std::thread::sleep(std::time::Duration::from_millis(500));
            }
            _ => {}
        }
        gfx.request_redraw();
    });
}
