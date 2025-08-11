mod gfx;
mod grid;
mod gridphx;
mod particles;
// mod sweep;
// mod sweepphx;
mod util;
mod vis;
use crate::vis::Visualization;
use std::thread;
fn main() {
    // let vis = Visualization::new();
    // vis.run();
    let a = thread::spawn(|| {
        let mut john = vec![0u32; 1000_000];
        for (i, j) in john.iter_mut().enumerate() {
            *j = i as u32;
        }
        dbg!(john[888]);
    });
    dbg!(&a);
    let a_result = a.join();
    dbg!(a_result);
}
