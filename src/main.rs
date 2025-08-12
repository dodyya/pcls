mod array;
mod gfx;
mod grid;
mod maybe_id;
mod particles;
mod phx;
mod vis;

use crate::vis::Visualization;
fn main() {
    let vis = Visualization::new();
    vis.run();
}
