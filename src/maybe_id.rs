use std::sync::atomic::{AtomicBool, AtomicUsize};

use crate::particles::O;

#[derive(Debug)]
pub struct MaybeID {
    some: AtomicBool,
    id: AtomicUsize,
}

impl MaybeID {
    pub fn is_some(&self) -> bool {
        self.some.load(O)
    }

    pub fn is_none(&self) -> bool {
        !self.some.load(O)
    }

    pub fn id(&self) -> Option<usize> {
        match self.some.load(O) {
            true => Some(self.id.load(O)),
            false => None,
        }
    }

    pub fn set(&self, val: usize) {
        self.some.store(true, O);
        self.id.store(val, O);
    }

    pub fn make_empty(&self) {
        self.some.store(false, O);
        self.id.store(0, O);
    }
}

impl Default for MaybeID {
    fn default() -> Self {
        MaybeID {
            some: AtomicBool::new(false),
            id: AtomicUsize::new(0),
        }
    }
}
