use crate::gfx::Gfx;
mod gfx;
fn main() {
    let mut g: Gfx = Gfx::new(600, 600);
    draw(&mut g.state, g.width.try_into().unwrap());
    g.run();
}
fn draw(state: &mut Vec<u8>, width: usize) {
    let start: Vec<i32> = vec![100, 200];
    let end: Vec<i32> = vec![200, 400];

    let mut x1 = start[0];
    let mut y1 = start[1];
    let x2 = end[0];
    let y2 = end[1];

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

        state[((y1 as usize * width as usize + x1 as usize) * 4)
            ..((y1 as usize * width as usize + x1 as usize) * 4 + 4)]
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
