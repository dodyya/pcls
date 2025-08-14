use crate::gfx;
use crate::sim::Simulation;
use atomic_float::AtomicF32;
use pixels::{Pixels, SurfaceTexture};
use rand::rngs::ThreadRng;
use rand::Rng;
use std::sync::atomic::Ordering::Relaxed as O;
use std::time::{Duration, Instant};
use winit::{
    dpi::PhysicalSize,
    event::{Event, MouseButton},
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

const MAX_PARTICLE_SIZE: f32 = 1.0 / 256.0;
const PARTICLES_ON_CLICK: usize = 250;
const WINDOW_SIZE: u32 = 1500;
const MASS: f32 = 1.0;

pub struct Visualization {
    window: Window,
    pixels: Pixels,
    sim: Simulation,
    event_loop: EventLoop<()>,
}
impl Visualization {
    pub fn run(mut self) {
        let mut cursor_pos: Option<(f32, f32)> = None;
        let mut last_frame = Instant::now();
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

            display(
                self.pixels.frame_mut(),
                self.sim.get_drawable(),
                WINDOW_SIZE,
            );
            _ = self.pixels.render();

            (frame_time, last_frame) = (last_frame.elapsed(), Instant::now());
            self.sim.step();

            if mouse_down {
                if let Some((cursor_x, cursor_y)) = cursor_pos {
                    if (0.0..=WINDOW_SIZE as f32).contains(&cursor_x)
                        && (0.0..=WINDOW_SIZE as f32).contains(&cursor_y)
                    {
                        add_particles(cursor_x, cursor_y, &mut self.sim, &mut rng);
                    }
                }
            }

            use winit::event::WindowEvent as we;
            match &event {
                Event::WindowEvent { event, .. } => match event {
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
                },
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
        let sim = Simulation::new(MAX_PARTICLE_SIZE * 2.0);

        Self {
            window,
            pixels,
            sim,
            event_loop,
        }
    }
}

fn display(frame: &mut [u8], particles: (&[AtomicF32], &[AtomicF32], &[AtomicF32]), width: u32) {
    gfx::_rst(frame);
    particles
        .0
        .iter()
        .zip(particles.1.iter())
        .zip(particles.2.iter())
        .for_each(|((x, y), r)| {
            gfx::draw_circle(frame, width, (x.load(O), y.load(O)), r.load(O));
        });
}

fn add_particles(cursor_x: f32, cursor_y: f32, sim: &mut Simulation, rng: &mut ThreadRng) {
    let sim_x = (cursor_x / WINDOW_SIZE as f32) * 2.0 - 1.0;
    let sim_y = 1.0 - (cursor_y / WINDOW_SIZE as f32) * 2.0;

    for _ in 0..PARTICLES_ON_CLICK {
        let dx = rng.gen_range(-0.2..0.2);
        let dy = rng.gen_range(-0.2..0.2);
        let r = MAX_PARTICLE_SIZE;

        sim.add_particle(sim_x + dx, sim_y + dy, r, 0.0, 0.0, MASS * r * r);
    }
}
