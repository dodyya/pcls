use pixels::{Pixels, SurfaceTexture};
use std::time::Instant;
use winit::{
    dpi::PhysicalSize,
    event::{Event, VirtualKeyCode, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

pub struct Gfx {
    pub window: Window,
    pixels: Pixels,
    width: u32,
    height: u32,
}

impl Gfx {
    pub fn new(width: u32, height: u32) -> (Self, EventLoop<()>) {
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

        (
            Gfx {
                window,
                pixels,
                width: size.width,
                height: size.height,
            },
            event_loop,
        )
    }

    pub fn draw_particles(&mut self, particles: (&[f32], &[f32], &[f32])) {
        particles
            .0
            .iter()
            .zip(particles.1.iter())
            .zip(particles.2.iter())
            .for_each(|((x, y), r)| {
                fill_circle(self.pixels.frame_mut(), self.width, (*x, *y), *r);
                draw_circle(self.pixels.frame_mut(), self.width, (*x, *y), *r);
            });
    }

    pub fn clear_frame(&mut self) {
        _rst(self.pixels.frame_mut());
    }

    pub fn render(&mut self) {
        self.pixels.render().unwrap();
    }

    pub fn request_redraw(&mut self) {
        self.window.request_redraw();
    }
}

fn _draw_line(frame: &mut [u8], width: u32, start: (i32, i32), end: (i32, i32)) {
    let (mut x1, mut y1) = start;
    let (x2, y2) = end;

    let dx = i32::abs(x2 - x1);
    let dy = i32::abs(y2 - y1);

    let sx = if x1 < x2 { 1 } else { -1 };
    let sy = if y1 < y2 { 1 } else { -1 };

    let mut err = dx - dy;
    let mut save_my_life = 10000;

    loop {
        if save_my_life < 0 {
            break;
        }

        _tpix(frame, width, (x1, y1));
        if x1 == x2 && y1 == y2 {
            break;
        }
        let e2 = 2 * err;
        if e2 > -dy {
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

pub fn draw_circle(frame: &mut [u8], width: u32, center: (f32, f32), radius: f32) {
    let pixel_center = ndc_to_pix(width, width, center);
    let pixel_radius = (radius * width as f32 / 2.0).round() as i32;
    _draw_circle(frame, width, pixel_center, pixel_radius);
}

fn _draw_circle(frame: &mut [u8], width: u32, center: (i32, i32), radius: i32) {
    let mut x = 0;
    let mut y = radius;
    let mut p = 1 - radius;

    while x <= y {
        if p < 0 {
            x += 1;
            p = p + 2 * x + 1;
        } else {
            x += 1;
            y -= 1;
            p = p + 2 * x - 2 * y + 1;
        }

        _tpix(frame, width, (x + center.0, y + center.1));
        _tpix(frame, width, (-x + center.0, y + center.1));
        _tpix(frame, width, (x + center.0, -y + center.1));
        _tpix(frame, width, (-x + center.0, -y + center.1));
        _tpix(frame, width, (y + center.0, x + center.1));
        _tpix(frame, width, (-y + center.0, x + center.1));
        _tpix(frame, width, (y + center.0, -x + center.1));
        _tpix(frame, width, (-y + center.0, -x + center.1));
    }
}

fn _fill_circle(frame: &mut [u8], width: u32, center: (i32, i32), radius: i32) {
    let mut x = 0;
    let mut y = radius;
    let mut p = 1 - radius;

    while x <= y {
        _crow(frame, width, (-x + center.0, x + center.0, y + center.1));
        _crow(frame, width, (-x + center.0, x + center.0, -y + center.1));
        _crow(frame, width, (-y + center.0, y + center.0, x + center.1));
        _crow(frame, width, (-y + center.0, y + center.0, -x + center.1));
        if p < 0 {
            x += 1;
            p = p + 2 * x + 1;
        } else {
            x += 1;
            y -= 1;
            p = p + 2 * x - 2 * y + 1;
        }
    }
}
pub fn fill_circle(frame: &mut [u8], width: u32, center: (f32, f32), radius: f32) {
    let pixel_center = ndc_to_pix(width, width, center);
    let pixel_radius = (radius * width as f32 / 2.0).round() as i32;
    _fill_circle(frame, width, pixel_center, pixel_radius);
}

pub fn draw_line(frame: &mut [u8], width: u32, height: u32, start: (f32, f32), end: (f32, f32)) {
    if f32::abs(start.0) > 1.0
        || f32::abs(start.1) > 1.0
        || f32::abs(end.0) > 1.0
        || f32::abs(end.1) > 1.0
    {
        return;
    }
    let pixel_start = ndc_to_pix(width, height, start);
    let pixel_end = ndc_to_pix(width, height, end);
    _draw_line(frame, width, pixel_start, pixel_end);
}

fn ndc_to_pix(width: u32, height: u32, (x, y): (f32, f32)) -> (i32, i32) {
    (
        (x * width as f32 / 2.0 + width as f32 / 2.0) as i32,
        (-y * height as f32 / 2.0 + height as f32 / 2.0) as i32,
    )
}
fn _tpix(frame: &mut [u8], width: u32, (x, y): (i32, i32)) {
    let index = (y.wrapping_mul(width as i32).wrapping_add(x) as usize).wrapping_mul(4);

    if index.saturating_add(4) <= frame.len() {
        let pixel = &mut frame[index..index + 4];
        if pixel == [0, 0, 0, 255] {
            pixel.copy_from_slice(&[255, 255, 255, 255]);
        } else {
            pixel.copy_from_slice(&[0, 0, 0, 255]);
        }
    } else {
        // println!("Tried to draw pixel at {:?}", (x, y))
    }
}

fn _cpix(frame: &mut [u8], width: u32, (x, y): (i32, i32)) {
    if x > width as i32 {
        println!("Tried to draw pixel at {:?}", (x, y));
    }
    let index = (y.wrapping_mul(width as i32).wrapping_add(x) as usize).wrapping_mul(4);

    if index.saturating_add(4) > frame.len() {
        println!("Tried to draw pixel at {:?}", (x, y));
        return;
    }
    let pixel = &mut frame[index..index + 4];
    pixel.copy_from_slice(&[255, 255, 255, 255]);
}

fn _crow(frame: &mut [u8], width: u32, (x1, x2, y): (i32, i32, i32)) {
    if x1 > x2
        || y < 0
        || x2 < 0
        || x1 >= width as i32
        || y >= (frame.len() / 4 / width as usize) as i32
    {
        return;
    }

    let x1 = x1.max(0);
    let x2 = x2.min(width as i32 - 1);

    let start_index = (y as usize * width as usize + x1 as usize) * 4;
    let end_index = (y as usize * width as usize + x2 as usize) * 4 + 4;
    let len = frame.len();

    if start_index >= len {
        return;
    }

    if end_index > len {
        let pixels = &mut frame[start_index..len];
        let white = vec![255; len - start_index];
        pixels.copy_from_slice(&white);
        return;
    }

    let pixels = &mut frame[start_index..end_index];
    let white = vec![255; end_index - start_index];
    pixels.copy_from_slice(&white);
}
fn _rst(frame: &mut [u8]) {
    // let black = [0, 0, 0, 255].repeat(frame.len() / 4);
    // frame.copy_from_slice(&black)
    frame.fill(0);
}
