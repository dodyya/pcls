mod gfx;
mod grid;
mod gridphx;
mod maybe_id;
mod particles;
// mod sweep;
// mod sweepphx;
mod array;
mod vis;
// use crate::gridphx::cell;
// use crate::util::Array3D;
use crate::vis::Visualization;
use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Mutex, RwLock};
use std::thread;
use std::time::Duration;
fn main() {
    let vis = Visualization::new();
    vis.run();
    // let state = RwLock::new(0u8);
    // dbg!(*state.read().unwrap());
    // *state.write().unwrap() = 5;
    // dbg!(*state.read().unwrap());
}
