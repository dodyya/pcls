use pixels::{Pixels, SurfaceTexture};
use winit::{
    dpi::PhysicalSize,
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

    pub fn run(mut self) {
        self.event_loop.run(move |event, _, control_flow| {
            println!("Event: {:?}", event);
            match &event {
                Event::WindowEvent { event, .. } => match event {
                    WindowEvent::CloseRequested => {
                        println!("Window closed!");
                        *control_flow = ControlFlow::Exit;
                    }
                    WindowEvent::CursorMoved { position, .. } => {
                        let x = position.x as usize;
                        let y = position.y as usize;

                        if x < self.width as usize && y < self.height as usize {
                            let index = (y * self.width as usize + x) * 4;
                            self.state[index..index + 4].copy_from_slice(&[0, 0, 0, 255]);
                        }
                    }
                    _ => {}
                },
                Event::RedrawRequested(_) => {
                    self.pixels.frame_mut().copy_from_slice(&self.state);
                    self.pixels.render().unwrap();
                }
                _ => {}
            }

            self.window.request_redraw();
        });
    }

    fn _draw_line(&mut self, start: (i32, i32), end: (i32, i32)) {
        let mut x1 = start.0;
        let mut y1 = start.1;
        let x2 = end.0;
        let y2 = end.1;

        let dx = i32::abs(x1 - x2);
        let dy = i32::abs(y1 - y2);

        let sx = if x1 < x2 { 1 } else { -1 };
        let sy = if y1 < y2 { 1 } else { -1 };

        let mut err = dx - dy;
        let mut save_my_life = 10000;

        loop {
            if save_my_life < 0 {
                break;
            }

            self.state[((y1 as usize * self.width as usize + x1 as usize) * 4)
                ..((y1 as usize * self.width as usize + x1 as usize) * 4 + 4)]
                .copy_from_slice(&[0, 0, 0, 255]);

            if x1 == x2 {
                break;
            }

            let e2 = 2 * err;
            if dy > -e2 {
                err -= dy;
                x1 += sx;
            }
            if e2 < dx {
                err += dx;
                y1 += sy;
            }

            save_my_life -= 1;
        }
    }
    pub fn draw_line(&mut self, start: (f32, f32), end: (f32, f32)) {
        if f32::abs(start.0) > 1.0
            || f32::abs(start.1) > 1.0
            || f32::abs(end.0) > 1.0
            || f32::abs(end.1) > 1.0
        {
            return;
        }
        let pixel_start = self.ndc_to_pixels(start);
        let pixel_end = self.ndc_to_pixels(end);
        self._draw_line(pixel_start, pixel_end);
    }
    fn ndc_to_pixels(&self, (x, y): (f32, f32)) -> (i32, i32) {
        (
            (x * self.width as f32 / 2.0 + self.width as f32 / 2.0) as i32,
            (-y * self.height as f32 / 2.0 + self.height as f32 / 2.0) as i32,
        )
    }
}
