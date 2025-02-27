use pixels::{Pixels, SurfaceTexture};
use winit::{
    dpi::LogicalSize,
    event::{Event, VirtualKeyCode, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

pub struct Gfx {
    event_loop: EventLoop<()>,
    window: Window,
    pixels: Pixels,
    /// The state; width\*height\*4 u8s representing RGBA
    pub state: Vec<u8>,
    /// Width of window in pixels
    pub width: u32,
    /// Height of window in pixels
    pub height: u32,
}

impl Gfx {
    /// Create a new Gfx instance
    pub fn new(width: u32, height: u32) -> Self {
        let event_loop = EventLoop::new();
        let window = WindowBuilder::new()
            .with_title("graphics")
            .with_inner_size(LogicalSize::new(width, height))
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

        let state = vec![255u8; size.width as usize * size.height as usize * 4];

        Gfx {
            event_loop,
            window,
            pixels,
            state,
            width: size.width,
            height: size.height,
        }
    }

    /// Run the event loop (Handles everything: input, rendering, updates)
    pub fn run(self) {
        let mut pixels = self.pixels; // Move pixels into the closure
        let state = self.state; // Move state into the closure
        let window = self.window; // Move window into the closure

        self.event_loop.run(move |event, _, control_flow| {
            match &event {
                Event::WindowEvent { event, .. } => match event {
                    WindowEvent::CloseRequested => {
                        println!("Window closed!");
                        *control_flow = ControlFlow::Exit;
                    }
                    WindowEvent::KeyboardInput { input, .. } => {
                        if let Some(VirtualKeyCode::Escape) = input.virtual_keycode {
                            println!("Escape key pressed! Exiting...");
                            *control_flow = ControlFlow::Exit;
                        }
                    }
                    _ => {}
                },
                Event::RedrawRequested(_) => {
                    pixels.frame_mut().copy_from_slice(&state);
                    pixels.render().unwrap();
                }
                _ => {}
            }

            window.request_redraw();
        });
    }
}
