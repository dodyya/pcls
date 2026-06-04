use crate::particles::XyHue;
use cuda_core::error::{DriverError, IntoResult};
use cuda_core::sys::{
    CUdeviceptr, CUgraphicsResource, cuGraphicsMapResources,
    cuGraphicsResourceGetMappedPointer_v2, cuGraphicsUnmapResources,
    cuGraphicsUnregisterResource,
};
use cuda_core::{CudaContext, CudaStream};
use glow::HasContext as _;
use std::ptr;
use std::sync::Arc;

// cuda.h (and so cuda-oxide's bindings) doesn't include cudaGL.h; one extern for the GL-
// specific entry point. The rest of GL interop is generic graphics-resource API already
// in cuda-oxide's `sys` re-export.
unsafe extern "C" {
    fn cuGraphicsGLRegisterBuffer(
        resource: *mut CUgraphicsResource,
        buffer: u32,
        flags: ::std::os::raw::c_uint,
    ) -> cuda_core::sys::CUresult;
}

// CU_GRAPHICS_REGISTER_FLAGS_WRITE_DISCARD: the extract kernel overwrites the whole buffer
// each frame; CUDA can skip reconciling prior contents.
const REGISTER_FLAGS_WRITE_DISCARD: u32 = 0x02;

// GL VBO holding `XyHue` triples, registered with CUDA so the sim's extract kernel writes
// into the same memory the renderer draws from. Map/unmap each frame: mapped → CUDA owns
// it, unmapped → GL owns it.
pub struct ParticleVbo {
    pub vbo: glow::Buffer,
    raw_name: u32,
    pub cap: usize,
    resource: CUgraphicsResource,
    mapped: Option<CUdeviceptr>,
    ctx: Arc<CudaContext>,
}

impl ParticleVbo {
    // Both the GL context and the CUDA context must be current on the calling thread.
    pub fn new(
        gl: &glow::Context,
        ctx: Arc<CudaContext>,
        cap: usize,
    ) -> Result<Self, DriverError> {
        ctx.bind_to_thread()?;
        // Safety: GL context is current per the function's contract.
        let (vbo, raw_name) = unsafe {
            let vbo = gl.create_buffer().expect("create vbo");
            let raw_name = vbo.0.get();
            gl.bind_buffer(glow::ARRAY_BUFFER, Some(vbo));
            gl.buffer_data_size(
                glow::ARRAY_BUFFER,
                (cap * size_of::<XyHue>()) as i32,
                glow::DYNAMIC_DRAW,
            );
            gl.bind_buffer(glow::ARRAY_BUFFER, None);
            (vbo, raw_name)
        };
        let mut resource: CUgraphicsResource = ptr::null_mut();
        // Safety: vbo is a valid GL buffer just created above; the GL context is current.
        unsafe {
            cuGraphicsGLRegisterBuffer(&mut resource, raw_name, REGISTER_FLAGS_WRITE_DISCARD)
                .result()?;
        }
        Ok(Self {
            vbo,
            raw_name,
            cap,
            resource,
            mapped: None,
            ctx,
        })
    }

    // Maps the VBO for CUDA access and returns the device pointer. The caller treats
    // it as a borrowed `[XyHue; cap]`; cuda-oxide's DisjointSlice (constructed at the
    // kernel-launch site) handles per-thread bounds checking.
    pub fn map(&mut self, stream: &CudaStream) -> Result<CUdeviceptr, DriverError> {
        debug_assert!(self.mapped.is_none(), "double map");
        // Safety: `self.resource` is registered; stream belongs to the same context.
        unsafe {
            cuGraphicsMapResources(1, &mut self.resource, stream.cu_stream()).result()?;
        }
        let mut devptr: CUdeviceptr = 0;
        let mut size: usize = 0;
        // Safety: the resource is currently mapped.
        unsafe {
            cuGraphicsResourceGetMappedPointer_v2(&mut devptr, &mut size, self.resource)
                .result()?;
        }
        self.mapped = Some(devptr);
        Ok(devptr)
    }

    pub fn unmap(&mut self, stream: &CudaStream) -> Result<(), DriverError> {
        debug_assert!(self.mapped.is_some(), "unmap without map");
        // Safety: `self.resource` was mapped on this stream.
        unsafe {
            cuGraphicsUnmapResources(1, &mut self.resource, stream.cu_stream()).result()?;
        }
        self.mapped = None;
        Ok(())
    }

    pub fn ensure_cap(
        &mut self,
        gl: &glow::Context,
        needed: usize,
    ) -> Result<(), DriverError> {
        if needed <= self.cap {
            return Ok(());
        }
        let new_cap = needed.next_power_of_two().max(64);
        debug_assert!(self.mapped.is_none(), "grow while mapped");
        self.ctx.bind_to_thread()?;
        // Safety: `self.resource` is registered; GL context is current.
        unsafe {
            cuGraphicsUnregisterResource(self.resource).result()?;
            gl.bind_buffer(glow::ARRAY_BUFFER, Some(self.vbo));
            gl.buffer_data_size(
                glow::ARRAY_BUFFER,
                (new_cap * size_of::<XyHue>()) as i32,
                glow::DYNAMIC_DRAW,
            );
            gl.bind_buffer(glow::ARRAY_BUFFER, None);
            cuGraphicsGLRegisterBuffer(
                &mut self.resource,
                self.raw_name,
                REGISTER_FLAGS_WRITE_DISCARD,
            )
            .result()?;
        }
        self.cap = new_cap;
        Ok(())
    }
}

impl Drop for ParticleVbo {
    fn drop(&mut self) {
        let _ = self.ctx.bind_to_thread();
        // Safety: `self.resource` is the resource we registered in new()/grow().
        unsafe {
            let _ = cuGraphicsUnregisterResource(self.resource);
        }
        // The GL buffer is freed by the caller's glow::Context teardown.
    }
}
