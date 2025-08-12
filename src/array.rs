#![allow(unused)]
use std::{
    fmt::{Debug, Display},
    ops::{Index, IndexMut},
    ptr::slice_from_raw_parts,
    sync::Arc,
    vec,
};

use crate::maybe_id::MaybeID;
#[derive(Clone, Debug)]
pub struct Array3D<T> {
    pub data: Vec<T>,
    pub width: usize,
    pub height: usize,
    pub depth: usize,
}

impl<T> Array3D<Option<T>> {
    pub fn new(width: usize, height: usize, depth: usize) -> Self {
        let mut data: Vec<Option<T>> = Vec::with_capacity(width * height * depth);
        for _ in 0..width * height * depth {
            data.push(None);
        }
        Array3D {
            data,
            width,
            height,
            depth,
        }
    }
}
impl<T> Array3D<T>
where
    T: Default + Copy,
{
    pub fn fill(data: T, width: usize, height: usize, depth: usize) -> Self {
        Array3D {
            data: vec![data; width * height * depth],
            width,
            height,
            depth,
        }
    }

    pub fn reset(&mut self, value: T) {
        self.data
            .copy_from_slice(&[value].repeat(self.width * self.height * self.depth));
    }
}

impl Array3D<MaybeID> {
    pub fn clear(&self) {
        for i in 0..self.data.len() {
            self.data[i].make_empty();
        }
    }
}

impl<T: Debug> Array3D<T> {
    pub fn from_vec(value: Vec<T>, width: usize, height: usize, depth: usize) -> Self {
        Array3D {
            data: value.try_into().unwrap(),
            width,
            height,
            depth,
        }
    }
}

impl<T> Array3D<T>
where
    T: Default,
{
    pub fn default(w: usize, h: usize, d: usize) -> Self {
        let mut data: Vec<T> = Vec::with_capacity(w * h * d);
        for i in 0..w * h * d {
            data.push(Default::default());
        }
        Array3D {
            data,
            width: w,
            height: h,
            depth: d,
        }
    }
}

impl<T> Index<(usize, usize, usize)> for Array3D<T> {
    type Output = T;

    fn index(&self, index: (usize, usize, usize)) -> &Self::Output {
        if index.0 >= self.width || index.1 >= self.height || index.2 >= self.depth {
            panic!(
                "Index out of bounds: ({},{},{}) not in [0,{})x[0,{})x[0,{})",
                index.0, index.1, index.2, self.width, self.height, self.depth
            );
        }
        &self.data[(index.0 * self.height + index.1) * self.depth + index.2]
    }
}

impl<T> IndexMut<(usize, usize, usize)> for Array3D<T> {
    fn index_mut(&mut self, index: (usize, usize, usize)) -> &mut Self::Output {
        if index.0 >= self.width || index.1 >= self.height || index.2 >= self.depth {
            panic!(
                "Index out of bounds: ({},{}, {}) not in [0,{})x[0,{})x[0,{}))",
                index.0, index.1, index.2, self.width, self.height, self.depth
            );
        }
        &mut self.data[(index.0 * self.height + index.1) * self.depth + index.2]
    }
}

impl<'a, T: 'a> Index<(usize, usize)> for Array3D<T> {
    type Output = [T];

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        if index.0 >= self.width || index.1 >= self.height {
            panic!(
                "Index out of bounds: ({},{}) not in [0,{})x[0,{})",
                index.0, index.1, self.width, self.height,
            );
        }
        let start = (index.0 * self.height + index.1) * self.depth;
        &self.data[start..start + self.depth]
    }
}

impl<'a, T: 'a> Index<(usize, usize)> for &Array3D<T> {
    type Output = [T];

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        if index.0 >= self.width || index.1 >= self.height {
            panic!(
                "Index out of bounds: ({},{}) not in [0,{})x[0,{})",
                index.0, index.1, self.width, self.height,
            );
        }
        let start = (index.0 * self.height + index.1) * self.depth;
        &self.data[start..start + self.depth]
    }
}

impl<'a, T: 'a> IndexMut<(usize, usize)> for Array3D<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        if index.0 >= self.width || index.1 >= self.height {
            panic!(
                "Index out of bounds: ({},{}) not in [0,{})x[0,{})",
                index.0, index.1, self.width, self.height,
            );
        }
        let start = (index.1 * self.width + index.0) * self.depth;
        &mut self.data[start..start + self.depth]
    }
}

impl<T> Display for Array3D<T>
where
    T: Display + Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for y in 0..self.height {
            writeln!(f)?;
            for x in 0..self.width {
                write!(f, "[")?;
                for d in 0..self.depth {
                    if d > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{}", Self::index(&self, (x, y, d)))?;
                }
                write!(f, "]")?;
            }
        }
        Ok(())
    }
}

impl<'a, T> Array3D<T> {
    pub unsafe fn split(&self, exp: u32) -> Vec<&[T]> {
        let n = 2usize.pow(exp);
        let mut out = Vec::with_capacity(n);
        let len = self.data.len() / n;
        let ptr = self.data.as_ptr();
        for i in 0..n {
            out.push(
                slice_from_raw_parts(ptr.add(i * len), len)
                    .as_ref()
                    .unwrap(),
            );
        }
        out
    }
}
