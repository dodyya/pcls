use crate::phx::Phx;
use crate::{gfx, particles::O};
use pixels::{Pixels, SurfaceTexture};
use std::time::{Duration, Instant};
use winit::{
    dpi::PhysicalSize,
    event::{Event, MouseButton},
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

use rand::Rng;
const MAX_PARTICLE_SIZE: f32 = 0.005;
const PARTICLE_COUNT: usize = 250;
const WINDOW_SIZE: u32 = 1500;
const MASS: f32 = 1.0;

pub struct Visualization {
    window: Window,
    pixels: Pixels,
    width: u32,
    height: u32,
    sim: Phx,
    event_loop: EventLoop<()>,
}
impl Visualization {
    pub fn run(mut self) {
        let mut cursor_pos: Option<(f32, f32)> = None;
        let mut last_frame_start = Instant::now();
        let mut frame_time = Duration::ZERO;
        let mut ticker: u8 = 0;
        let mut rng = rand::thread_rng();
        let mut mouse_down = false;

        self.event_loop.run(move |event, _, control_flow| {
            if ticker == 0 {
                self.window.set_title(&format!(
                    "Verlet particle simulation: {} particles - FPS: {:.0}",
                    self.sim.pcls.count,
                    1.0 / frame_time.as_secs_f64() as f64
                ));
            }
            ticker = ticker.wrapping_add(8);

            display(self.pixels.frame_mut(), &self.sim, self.width);
            _ = self.pixels.render();

            frame_time = last_frame_start.elapsed();
            last_frame_start = Instant::now();
            self.sim.step();

            if mouse_down {
                if let Some((cursor_x, cursor_y)) = cursor_pos {
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
                            let r = MAX_PARTICLE_SIZE;

                            self.sim.add_particle(
                                sim_x + dx,
                                sim_y + dy,
                                r,
                                0.0,
                                0.0,
                                MASS * r * r,
                            );
                            if PARTICLE_COUNT == 1 {
                                mouse_down = false;
                            }
                        }
                    }
                }
            }

            use winit::event::WindowEvent as we;
            match &event {
                Event::WindowEvent { window_id, event } => {
                    if *window_id != self.window.id() {
                        return;
                    }
                    match event {
                        we::CloseRequested => {
                            *control_flow = ControlFlow::Exit;
                            return;
                        }
                        we::CursorMoved { position, .. } => {
                            cursor_pos = Some((position.x as f32, position.y as f32));
                        }
                        we::MouseInput {
                            state: winit::event::ElementState::Pressed,
                            button,
                            ..
                        } => {
                            if *button == MouseButton::Left {
                                mouse_down = true;
                            } else if *button == MouseButton::Right {
                                self.sim.clear();
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
                                    virtual_keycode: Some(k),
                                    state: winit::event::ElementState::Pressed,
                                    ..
                                },
                            ..
                        } => match k {
                            winit::event::VirtualKeyCode::Space => {
                                self.sim.toggle_gravity();
                            }
                            winit::event::VirtualKeyCode::S => {
                                self.sim.stop();
                            }
                            _ => {}
                        },
                        _ => {}
                    }
                }
                _ => {}
            }
        });
    }

    pub fn new() -> Self {
        let width = WINDOW_SIZE;
        let height = WINDOW_SIZE;
        let event_loop = EventLoop::new();
        let window = WindowBuilder::new()
            .with_title("graphics")
            .with_inner_size(PhysicalSize::new(width, height))
            .with_resizable(false)
            .build(&event_loop)
            .unwrap();

        let size = window.inner_size();
        let pixels = Pixels::new(
            size.width,
            size.height,
            SurfaceTexture::new(size.width, size.height, &window),
        )
        .unwrap();

        Self {
            window,
            pixels,
            width,
            height,
            sim: Phx::new(MAX_PARTICLE_SIZE * 2.0),
            event_loop,
        }
    }
}
fn display(frame: &mut [u8], sim: &Phx, width: u32) {
    gfx::_rst(frame);
    let particles = sim.get_drawable_particles();
    particles
        .0
        .iter()
        .zip(particles.1.iter())
        .zip(particles.2.iter())
        .for_each(|((x, y), r)| {
            gfx::fill_circle(frame, width, (x.load(O), y.load(O)), r.load(O));
            gfx::draw_circle(frame, width, (x.load(O), y.load(O)), r.load(O));
        });
}
