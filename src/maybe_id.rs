use std::sync::atomic::Ordering::Relaxed as O;
use std::sync::atomic::{AtomicBool, AtomicUsize};

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

    pub fn unchecked_id(&self) -> usize {
        self.id.load(O)
    }

    pub fn set(&self, val: usize) {
        self.some.store(true, O);
        self.id.store(val, O);
    }

    pub fn make_empty(&self) {
        self.some.store(false, O);
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
