mod array;
mod gfx;
mod grid;
mod maybe_id;
mod particles;
mod sim;
mod vis;

use crate::vis::Visualization;
fn main() {
    let vis = Visualization::new();
    vis.run();
}
