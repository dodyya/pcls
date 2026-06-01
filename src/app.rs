use crate::render::Renderer;
use crate::sim::Simulation;
use crate::types::Real;
use glutin::{
    config::ConfigTemplateBuilder,
    context::{ContextAttributesBuilder, PossiblyCurrentContext},
    display::GetGlDisplay,
    prelude::*,
    surface::{Surface, SurfaceAttributesBuilder, WindowSurface},
};
use glutin_winit::{DisplayBuilder, GlWindow};
use rand::Rng;
use rand::rngs::ThreadRng;
use raw_window_handle::HasWindowHandle;
use std::num::NonZeroU32;
use std::time::{Duration, Instant};
use winit::{
    application::ApplicationHandler,
    dpi::PhysicalSize,
    event::MouseButton,
    event_loop::{ActiveEventLoop, ControlFlow, EventLoop},
    keyboard::PhysicalKey,
    window::Window,
};

const MAX_PARTICLE_SIZE: Real = 1. / 256.;
const PARTICLES_ON_CLICK: usize = 250;
const WINDOW_SIZE: u32 = 1500;
const PDENSITY: Real = 1.0;

pub struct App {
    vis: Option<Visualization>,
}

impl App {
    pub fn run() {
        let event_loop = EventLoop::new().unwrap();
        event_loop.set_control_flow(ControlFlow::Poll);
        let mut app = Self { vis: None };
        let _ = event_loop.run_app(&mut app);
    }
}

struct Visualization {
    sim: Simulation,
    st: SimState,
    gl_context: PossiblyCurrentContext,
    gl_surface: Surface<WindowSurface>,
    renderer: Renderer,
    window: Window,
}

struct SimState {
    cursor_pos: Option<(Real, Real)>,
    last_frame: Instant,
    frame_time: Duration,
    ticker: u8,
    rng: ThreadRng,
    mouse_down: bool,
}

impl SimState {
    fn new() -> Self {
        Self {
            cursor_pos: None,
            last_frame: Instant::now(),
            frame_time: Duration::ZERO,
            ticker: 0,
            rng: rand::thread_rng(),
            mouse_down: false,
        }
    }
}

impl Visualization {
    fn new(event_loop: &ActiveEventLoop) -> Self {
        let window_attributes = Window::default_attributes()
            .with_title("graphics")
            .with_inner_size(PhysicalSize::new(WINDOW_SIZE, WINDOW_SIZE))
            .with_resizable(false);

        // glutin-winit handles the platform-specific window-vs-config ordering.
        let display_builder =
            DisplayBuilder::new().with_window_attributes(Some(window_attributes));
        let (window, gl_config) = display_builder
            .build(event_loop, ConfigTemplateBuilder::new(), |configs| {
                configs
                    .reduce(|acc, cfg| {
                        if cfg.num_samples() > acc.num_samples() {
                            cfg
                        } else {
                            acc
                        }
                    })
                    .unwrap()
            })
            .unwrap();
        let window = window.unwrap();

        let raw_window_handle = window.window_handle().unwrap().as_raw();
        let gl_display = gl_config.display();

        let context_attributes =
            ContextAttributesBuilder::new().build(Some(raw_window_handle));
        // SAFETY: `raw_window_handle` is valid because `window` is kept alive in `Visualization`
        // for the lifetime of the context (struct field order ensures context drops before window).
        let not_current_context = unsafe {
            gl_display
                .create_context(&gl_config, &context_attributes)
                .unwrap()
        };

        let surface_attributes =
            window.build_surface_attributes(SurfaceAttributesBuilder::new()).unwrap();
        // SAFETY: same as above — `window` outlives the surface.
        let gl_surface = unsafe {
            gl_display
                .create_window_surface(&gl_config, &surface_attributes)
                .unwrap()
        };

        let gl_context = not_current_context.make_current(&gl_surface).unwrap();

        // SAFETY: the context was made current just above, so `get_proc_address` returns
        // valid GL function pointers compatible with the calls glow will make through them.
        let gl = unsafe {
            glow::Context::from_loader_function_cstr(|s| gl_display.get_proc_address(s).cast())
        };
        let renderer = Renderer::new(gl);

        Self {
            sim: Simulation::new(MAX_PARTICLE_SIZE * 2.0),
            st: SimState::new(),
            gl_context,
            gl_surface,
            renderer,
            window,
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.vis.is_none() {
            self.vis = Some(Visualization::new(event_loop));
        }
    }

    fn window_event(
        &mut self,
        event_loop: &ActiveEventLoop,
        _window_id: winit::window::WindowId,
        event: winit::event::WindowEvent,
    ) {
        use winit::event::WindowEvent as we;
        use winit::keyboard::KeyCode as kc;

        let Some(vis) = self.vis.as_mut() else {
            return;
        };

        match event {
            we::CloseRequested => {
                event_loop.exit();
            }
            we::Resized(size) => {
                if let (Some(w), Some(h)) =
                    (NonZeroU32::new(size.width), NonZeroU32::new(size.height))
                {
                    vis.gl_surface.resize(&vis.gl_context, w, h);
                }
            }
            we::RedrawRequested => {
                vis.st.frame_time = vis.st.last_frame.elapsed();
                vis.st.last_frame = Instant::now();

                if vis.st.mouse_down {
                    if let Some((cursor_x, cursor_y)) = vis.st.cursor_pos {
                        if (0.0..=WINDOW_SIZE as Real).contains(&cursor_x)
                            && (0.0..=WINDOW_SIZE as Real).contains(&cursor_y)
                        {
                            add_particles(cursor_x, cursor_y, &mut vis.sim, &mut vis.st.rng);
                        }
                    }
                }

                if vis.st.ticker % 16 == 0 {
                    vis.window.set_title(&format!(
                        "Verlet particle simulation: {} particles - FPS: {:.0}",
                        vis.sim.pcls.count,
                        1.0 / vis.st.frame_time.as_secs_f64(),
                    ));
                }
                vis.st.ticker = vis.st.ticker.wrapping_add(1);

                vis.sim.step();

                vis.renderer.render(
                    vis.sim.get_drawable(),
                    vis.sim.is_hue_force_enabled(),
                    vis.window.inner_size().width,
                );
                vis.gl_surface.swap_buffers(&vis.gl_context).unwrap();
            }
            we::CursorMoved { position, .. } => {
                vis.st.cursor_pos = Some((position.x as Real, position.y as Real));
            }
            we::MouseInput {
                state: winit::event::ElementState::Pressed,
                button,
                ..
            } => {
                if button == MouseButton::Left {
                    vis.st.mouse_down = true;
                } else if button == MouseButton::Right {
                    vis.sim.clear();
                }
            }
            we::MouseInput {
                state: winit::event::ElementState::Released,
                button: winit::event::MouseButton::Left,
                ..
            } => {
                vis.st.mouse_down = false;
            }
            we::KeyboardInput {
                event:
                    winit::event::KeyEvent {
                        physical_key: PhysicalKey::Code(k),
                        state: winit::event::ElementState::Pressed,
                        ..
                    },
                ..
            } => match k {
                kc::Space => {
                    vis.sim.toggle_gravity();
                }
                kc::KeyM => {
                    vis.sim.toggle_hue_force();
                }
                kc::KeyD => {
                    vis.sim.toggle_donut();
                }
                kc::KeyS => {
                    vis.sim.stop();
                }
                kc::KeyQ => event_loop.exit(),
                _ => {}
            },
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _event_loop: &ActiveEventLoop) {
        if let Some(vis) = &self.vis {
            vis.window.request_redraw();
        }
    }

    fn exiting(&mut self, _event_loop: &ActiveEventLoop) {
        // Drop GL surface/context (and the window) while the display is still alive.
        self.vis = None;
    }
}

fn add_particles(cursor_x: Real, cursor_y: Real, sim: &mut Simulation, rng: &mut ThreadRng) {
    let sim_x = (cursor_x / WINDOW_SIZE as Real) * 2.0 - 1.0;
    let sim_y = 1.0 - (cursor_y / WINDOW_SIZE as Real) * 2.0;

    for _ in 0..PARTICLES_ON_CLICK {
        let dx: Real = rng.gen_range(-0.2..0.2);
        let dy: Real = rng.gen_range(-0.2..0.2);
        let r = MAX_PARTICLE_SIZE;
        let hue: Real = rng.gen_range(0.0..1.0);
        sim.add_particle(sim_x + dx, sim_y + dy, r, PDENSITY * r * r, hue);
    }
}
