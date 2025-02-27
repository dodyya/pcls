use crate::gfx::Gfx;
mod gfx;
fn main() {
    let mut g: Gfx = Gfx::new(600, 600);
    draw(&mut g.state, g.width.try_into().unwrap());
    g.run();
}
fn draw(state: &mut Vec<u8>, width: usize) {
    for i in 0..width {
        let index = (i * width + i) * 4;
        state[index..index + 4].copy_from_slice(&[0, 0, 0, 255]);
    }
}
