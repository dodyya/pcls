use crate::gfx::Gfx;
mod gfx;
fn main() {
    let mut g: Gfx = Gfx::new(600, 600);
    g.draw_line((-0.65, -0.2), (0.88, 0.7));
    g.run();
}
