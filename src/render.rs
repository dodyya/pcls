use crate::gl_interop::ParticleVbo;
use crate::sim::PARTICLE_RADIUS;
use crate::types::Real;
use glow::HasContext as _;

const VERTEX_SHADER: &str = r#"#version 330 core
layout (location = 0) in vec2 a_pos;
layout (location = 1) in float a_hue;
uniform float u_point_size;
out float v_hue;
void main() {
    gl_Position = vec4(a_pos, 0.0, 1.0);
    gl_PointSize = u_point_size;
    v_hue = a_hue;
}
"#;

// HSV->RGB by Iñigo Quilez (saturation = value = 1, hue periodic via fract).
const FRAGMENT_SHADER: &str = r#"#version 330 core
in float v_hue;
uniform int u_hue_force;
out vec4 frag_color;
vec3 hue_to_rgb(float h) {
    vec3 k = vec3(1.0, 2.0 / 3.0, 1.0 / 3.0);
    vec3 p = abs(fract(vec3(h) + k) * 6.0 - vec3(3.0));
    return clamp(p - vec3(1.0), 0.0, 1.0);
}
void main() {
    vec2 d = gl_PointCoord - vec2(0.5);
    if (dot(d, d) > 0.25) {
        discard;
    }
    vec3 color = u_hue_force != 0 ? hue_to_rgb(v_hue) : vec3(1.0);
    frag_color = vec4(color, 1.0);
}
"#;

const FLOATS_PER_PARTICLE: usize = 3; // x, y, hue

pub struct Renderer {
    program: glow::Program,
    vao: glow::VertexArray,
    u_point_size: Option<glow::UniformLocation>,
    u_hue_force: Option<glow::UniformLocation>,
}

impl Renderer {
    // Binds the VAO's vertex attribute pointers to `vbo`. Caller guarantees the GL context
    // is current; `vbo`'s GL name must outlive this Renderer.
    pub fn new(gl: &glow::Context, vbo: &ParticleVbo) -> Self {
        // Safety: caller guarantees a GL context is current on this thread.
        unsafe {
            let program = link_program(gl, VERTEX_SHADER, FRAGMENT_SHADER);

            let vao = gl.create_vertex_array().unwrap();
            gl.bind_vertex_array(Some(vao));
            gl.bind_buffer(glow::ARRAY_BUFFER, Some(vbo.vbo));

            let stride = (FLOATS_PER_PARTICLE * size_of::<f32>()) as i32;
            gl.enable_vertex_attrib_array(0);
            gl.vertex_attrib_pointer_f32(0, 2, glow::FLOAT, false, stride, 0);
            gl.enable_vertex_attrib_array(1);
            gl.vertex_attrib_pointer_f32(
                1,
                1,
                glow::FLOAT,
                false,
                stride,
                2 * size_of::<f32>() as i32,
            );

            gl.bind_vertex_array(None);
            gl.bind_buffer(glow::ARRAY_BUFFER, None);

            let u_point_size = gl.get_uniform_location(program, "u_point_size");
            let u_hue_force = gl.get_uniform_location(program, "u_hue_force");

            gl.enable(glow::PROGRAM_POINT_SIZE);
            gl.clear_color(0.0, 0.0, 0.0, 1.0);

            Self {
                program,
                vao,
                u_point_size,
                u_hue_force,
            }
        }
    }

    pub fn render(
        &self,
        gl: &glow::Context,
        count: usize,
        hue_force_enabled: bool,
        win_size: u32,
    ) {
        let point_size = (PARTICLE_RADIUS * win_size as Real) as _;
        // Safety: GL context is current; VAO/VBO/program are live.
        unsafe {
            gl.viewport(0, 0, win_size as i32, win_size as i32);
            gl.clear(glow::COLOR_BUFFER_BIT);
            if count > 0 {
                gl.use_program(Some(self.program));
                gl.uniform_1_f32(self.u_point_size.as_ref(), point_size);
                gl.uniform_1_i32(self.u_hue_force.as_ref(), hue_force_enabled as i32);
                gl.bind_vertex_array(Some(self.vao));
                gl.draw_arrays(glow::POINTS, 0, count as i32);
                gl.bind_vertex_array(None);
            }
        }
    }
}

/// # Safety
/// A GL context must be current on the calling thread.
unsafe fn link_program(gl: &glow::Context, vs_src: &str, fs_src: &str) -> glow::Program {
    // Safety: GL context is current per this function's contract.
    unsafe {
        let program = gl.create_program().unwrap();
        let shaders = [(glow::VERTEX_SHADER, vs_src), (glow::FRAGMENT_SHADER, fs_src)];
        let mut compiled = Vec::with_capacity(shaders.len());
        for (kind, src) in shaders {
            let shader = gl.create_shader(kind).unwrap();
            gl.shader_source(shader, src);
            gl.compile_shader(shader);
            if !gl.get_shader_compile_status(shader) {
                panic!("shader compile error: {}", gl.get_shader_info_log(shader));
            }
            gl.attach_shader(program, shader);
            compiled.push(shader);
        }
        gl.link_program(program);
        if !gl.get_program_link_status(program) {
            panic!("program link error: {}", gl.get_program_info_log(program));
        }
        for shader in compiled {
            gl.detach_shader(program, shader);
            gl.delete_shader(shader);
        }
        program
    }
}
