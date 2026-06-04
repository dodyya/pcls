use cuda_core::DeviceBuffer;
use std::marker::PhantomData;
use std::mem::ManuallyDrop;
use std::ops::{Deref, DerefMut};

// Non-owning, length-bounded view over a DeviceBuffer. Wraps a ManuallyDrop<DeviceBuffer>
// so its Drop never calls cuMemFree, and deref-coerces to &DeviceBuffer<T> / &mut
// DeviceBuffer<T> for the cuda-host sync launch macros. Replaces the raw
// `ManuallyDrop::new(DeviceBuffer::from_raw_parts(...))` pattern with a borrow-typed
// handle that cannot accidentally free the underlying allocation.
//
// The Arc clone of the source buffer's context still happens at construction time —
// cuda-oxide's sync launch path requires a `&DeviceBuffer<T>`, which carries an
// Arc<CudaContext>. Removing that cost is an upstream change.
pub struct DeviceView<'a, T> {
    buf: ManuallyDrop<DeviceBuffer<T>>,
    _p: PhantomData<&'a mut DeviceBuffer<T>>,
}

impl<'a, T> DeviceView<'a, T> {
    pub fn from_buffer(buf: &'a DeviceBuffer<T>, len: usize) -> Self {
        assert!(len <= buf.len(), "DeviceView len exceeds buffer len");
        // Safety: ManuallyDrop suppresses the inner DeviceBuffer's cuMemFree; the source
        // `buf` (lifetime-bound by `'a`) keeps ownership.
        let inner =
            unsafe { DeviceBuffer::from_raw_parts(buf.cu_deviceptr(), len, buf.context().clone()) };
        Self {
            buf: ManuallyDrop::new(inner),
            _p: PhantomData,
        }
    }
}

impl<'a, T> Deref for DeviceView<'a, T> {
    type Target = DeviceBuffer<T>;
    fn deref(&self) -> &DeviceBuffer<T> {
        &self.buf
    }
}

impl<'a, T> DerefMut for DeviceView<'a, T> {
    fn deref_mut(&mut self) -> &mut DeviceBuffer<T> {
        &mut self.buf
    }
}
