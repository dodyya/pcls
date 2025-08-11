mod gfx;
mod grid;
mod gridphx;
mod particles;
// mod sweep;
// mod sweepphx;
mod util;
mod vis;
use crate::gridphx::cell;
use crate::util::Array3D;
use crate::vis::Visualization;
use std::cell::RefCell;
use std::rc::Rc;
use std::sync::{Arc, Mutex, RwLock};
use std::thread;
fn main() {
    // let vis = Visualization::new();
    // vis.run();
    let mut arr = Array3D::<usize>::new(2, 3, 4);
    for i in 0..2 {
        for j in 0..3 {
            for k in 0..4 {
                arr[(i, j, k)] = (i + 1) * 100 + (j + 1) * 10 + (k + 1);
            }
        }
    }
    //Should display 3 rows, 2 cols, each cell 4 deep
    //[111, 112, 113, 114] [211, 212, 213, 214]
    //[121, 122, 123, 124] [221, 222, 223, 224]
    //[131, 132, 133, 134] [231, 232, 233, 234]
    // println!("{}", arr);
    let state: Arc<Vec<RwLock<u32>>> = Arc::new((0..300).map(|v| RwLock::new(v)).collect());
    unsafe {
        let split = arr.split(1);
        let a = split[0];
        let b = split[1];
        let state_a = Arc::clone(&state);
        let state_b = Arc::clone(&state);
        let scope = thread::scope(move |s| {
            let a_thread = s.spawn(move || {
                for i in 0..1 {
                    for j in 0..3 {
                        for k in 0..4 {
                            let data = a[ind(i, j, k, 3, 4)];
                            *state_a[data].write().unwrap() = 1;
                        }
                    }
                }
            });
            let b_thread = s.spawn(move || {
                for i in 0..1 {
                    for j in 0..3 {
                        for k in 0..4 {
                            let data = b[ind(i, j, k, 3, 4)];
                            *state_b[data].write().unwrap() = 2;
                        }
                    }
                }
            });
        });

        dbg!(state);
    }

    fn ind(i: usize, j: usize, k: usize, height: usize, depth: usize) -> usize {
        (i * height + j) * depth + k
    }
}
