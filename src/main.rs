use std::time::{Duration, Instant};

use winit::event::WindowEvent as we;
use winit::{
    event::{ElementState, Event, MouseButton},
    event_loop::ControlFlow,
};

mod gfx;
mod grid;
mod gridphx;
mod particles;
mod sweep;
mod sweepphx;
use gfx::Gfx;
use gridphx::Phx;
use rand::Rng;
// use sweepphx::Phx;
const MAX_PARTICLE_SIZE: f32 = 0.01;
const PARTICLE_COUNT: usize = 25;
const WINDOW_SIZE: u32 = 1500;
const MASS: f32 = 1.0;

fn main() {
    let (mut gfx, event_loop) = Gfx::new(WINDOW_SIZE, WINDOW_SIZE);

    let mut sim = Phx::new(MAX_PARTICLE_SIZE * 2.0);

    let mut last_cursor_pos: Option<(f32, f32)> = None;
    let mut mouse_down = false;
    let mut rng = rand::thread_rng();
    let mut ticker: u8 = 0;
    let mut last_frame_start = Instant::now();
    let mut frame_time = Duration::ZERO;

    event_loop.run(move |event, _, control_flow| {
        if ticker == 0 {
            gfx.window.set_title(&format!(
                "Verlet particle simulation: {} particles - FPS: {:.0}",
                sim.pcls.count,
                1.0 / frame_time.as_secs_f64() as f64
            ))
        }
        ticker = ticker.wrapping_add(8);

        gfx.clear_frame();
        gfx.draw_particles(sim.get_drawable_particles());
        gfx.render();

        frame_time = last_frame_start.elapsed();
        last_frame_start = Instant::now();
        sim.step();

        if mouse_down {
            if let Some((cursor_x, cursor_y)) = last_cursor_pos {
                if cursor_x >= 0.0
                    && cursor_x <= WINDOW_SIZE as f32
                    && cursor_y >= 0.0
                    && cursor_y <= WINDOW_SIZE as f32
                {
                    let sim_x = (cursor_x / WINDOW_SIZE as f32) * 2.0 - 1.0;
                    let sim_y = 1.0 - (cursor_y / WINDOW_SIZE as f32) * 2.0;

                    for _ in 0..PARTICLE_COUNT {
                        let dx = rng.gen_range(-0.2..0.2);
                        let dy = rng.gen_range(-0.2..0.2);
                        let r = rng.gen_range(0.005..MAX_PARTICLE_SIZE);

                        sim.add_particle(sim_x + dx, sim_y + dy, r, 0.0, 0.0, MASS * r * r);
                    }
                }
            }
        }

        match &event {
            Event::WindowEvent { window_id, event } => {
                if *window_id != gfx.window.id() {
                    return;
                }
                match event {
                    we::CloseRequested => {
                        *control_flow = ControlFlow::Exit;
                        return;
                    }
                    we::CursorMoved { position, .. } => {
                        last_cursor_pos = Some((position.x as f32, position.y as f32));
                    }
                    we::MouseInput {
                        state: winit::event::ElementState::Pressed,
                        button,
                        ..
                    } => {
                        if *button == MouseButton::Left {
                            mouse_down = true;
                        } else if *button == MouseButton::Right {
                            sim.clear();
                        }
                    }
                    we::MouseInput {
                        state: winit::event::ElementState::Released,
                        button: winit::event::MouseButton::Left,
                        ..
                    } => {
                        mouse_down = false;
                    }
                    we::KeyboardInput {
                        input:
                            winit::event::KeyboardInput {
                                virtual_keycode: Some(winit::event::VirtualKeyCode::Space),
                                state: winit::event::ElementState::Pressed,
                                ..
                            },
                        ..
                    } => {
                        sim.toggle_gravity();
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    });
}
