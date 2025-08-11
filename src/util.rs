#![allow(unused)]
use std::{
    fmt::{Debug, Display},
    ops::{Index, IndexMut},
    vec,
};

#[derive(Clone)]
pub struct Array2D<T> {
    pub data: Vec<T>,
    pub width: usize,
    pub height: usize,
}

impl<T> Array2D<T>
where
    T: Clone + Default + Copy,
{
    pub fn new(width: usize, height: usize) -> Self {
        Array2D {
            data: vec![Default::default(); width * height],
            width,
            height,
        }
    }

    pub fn fill(data: T, width: usize, height: usize) -> Self {
        Array2D {
            data: vec![data; width * height],
            width,
            height,
        }
    }

    pub fn zero(&mut self) {
        self.data
            .copy_from_slice(&[Default::default()].repeat(self.width * self.height));
    }

    pub fn reset(&mut self, value: T) {
        self.data
            .copy_from_slice(&[value].repeat(self.width * self.height));
    }
}

impl<T: Debug> Array2D<T> {
    pub fn from_vec(value: Vec<T>, width: usize, height: usize) -> Self {
        Array2D {
            data: value.try_into().unwrap(),
            width,
            height,
        }
    }
}

impl<T> Array2D<T>
where
    T: Default,
{
    pub fn default(w: usize, h: usize) -> Self {
        let mut data: Vec<T> = Vec::with_capacity(w * h);
        for _ in 0..w * h {
            data.push(Default::default());
        }
        Array2D {
            data,
            width: w,
            height: h,
        }
    }
}

impl<T> Index<(usize, usize)> for Array2D<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        if index.0 >= self.width || index.1 >= self.height {
            panic!(
                "Index out of bounds: ({},{}) not in [0,{})x[0,{}))",
                index.0, index.1, self.width, self.height
            );
        }
        &self.data[index.1 * self.width + index.0]
    }
}

impl<T> IndexMut<(usize, usize)> for Array2D<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        if index.0 >= self.width || index.1 >= self.height {
            panic!(
                "Index out of bounds: ({},{}) not in [0,{})x[0,{}))",
                index.0, index.1, self.width, self.height
            );
        }
        &mut self.data[index.1 * self.width + index.0]
    }
}

impl<T> Display for Array2D<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for y in 0..self.height {
            writeln!(f)?;
            for x in 0..self.width {
                write!(f, "{:^9}", format!("{} ", self[(x, y)]))?;
            }
        }
        Ok(())
    }
}

impl<T> Debug for Array2D<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self)?;
        Ok(())
    }
}

impl<T> Array2D<T>
where
    T: Display + Copy,
{
    pub(crate) fn fill_circle(&mut self, center_x: i32, center_y: i32, radius: f32, value: T) {
        let r_squared = (radius * radius) as i32;

        for y in (center_y - radius.ceil() as i32)..=(center_y + radius.ceil() as i32) {
            for x in (center_x - radius.ceil() as i32)..=(center_x + radius.ceil() as i32) {
                if x >= 0 && y >= 0 && x < self.width as i32 && y < self.height as i32 {
                    let dx = x - center_x;
                    let dy = y - center_y;

                    if dx * dx + dy * dy <= r_squared {
                        self[(x as usize, y as usize)] = value;
                    }
                }
            }
        }
    }
}

// ----
#[derive(Clone, Debug)]
pub struct Array3D<T> {
    pub data: Vec<T>,
    pub width: usize,
    pub height: usize,
    pub depth: usize,
}

impl<T> Array3D<T>
where
    T: Clone + Default + Copy,
{
    pub fn new(width: usize, height: usize, depth: usize) -> Self {
        Array3D {
            data: vec![Default::default(); width * height * depth],
            width,
            height,
            depth,
        }
    }

    pub fn fill(data: T, width: usize, height: usize, depth: usize) -> Self {
        Array3D {
            data: vec![data; width * height * depth],
            width,
            height,
            depth,
        }
    }

    pub fn clear(&mut self) {
        self.data.fill(Default::default());
    }

    pub fn reset(&mut self, value: T) {
        self.data
            .copy_from_slice(&[value].repeat(self.width * self.height * self.depth));
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
        &self.data[(index.1 * self.width + index.0) * self.height + index.2]
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
        &mut self.data[(index.1 * self.width + index.0) * self.height + index.2]
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
        let start = (index.1 * self.width + index.0) * self.depth;
        &self.data[start..start + self.depth]
    }
}

impl<'a, T: 'a> IndexMut<(usize, usize)> for Array3D<T> {
    // type Output = [T];

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
                write!(f, "{:^9}", format!("{:?} ", &self[(x, y)]))?;
            }
        }
        Ok(())
    }
}

// impl<T> Debug for Array3D<T>
// where
//     T: Display + Debug,
// {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         write!(f, "{}", self)?;
//         Ok(())
//     }
// }
