use glow::HasContext as _;

const VERTEX_SHADER: &str = r#"#version 330 core
layout (location = 0) in vec2 a_pos;
layout (location = 1) in float a_radius;
layout (location = 2) in float a_charge;
uniform float u_win_size;
out float v_charge;
void main() {
    gl_Position = vec4(a_pos, 0.0, 1.0);
    // NDC radius r spans r * (win/2) px; diameter (point size) = r * win.
    gl_PointSize = a_radius * u_win_size;
    v_charge = a_charge;
}
"#;

const FRAGMENT_SHADER: &str = r#"#version 330 core
in float v_charge;
uniform int u_coulomb;
out vec4 frag_color;
void main() {
    vec2 d = gl_PointCoord - vec2(0.5);
    if (dot(d, d) > 0.25) {
        discard;
    }
    vec3 color;
    if (u_coulomb != 0) {
        color = v_charge > 0.0 ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 0.0, 1.0);
    } else {
        color = vec3(1.0);
    }
    frag_color = vec4(color, 1.0);
}
"#;

const FLOATS_PER_PARTICLE: usize = 4; // x, y, radius, charge

pub struct Renderer {
    gl: glow::Context,
    program: glow::Program,
    vao: glow::VertexArray,
    vbo: glow::Buffer,
    u_win_size: Option<glow::UniformLocation>,
    u_coulomb: Option<glow::UniformLocation>,
    verts: Vec<f32>,
}

impl Renderer {
    pub fn new(gl: glow::Context) -> Self {
        unsafe {
            let program = link_program(&gl, VERTEX_SHADER, FRAGMENT_SHADER);

            let vao = gl.create_vertex_array().unwrap();
            let vbo = gl.create_buffer().unwrap();
            gl.bind_vertex_array(Some(vao));
            gl.bind_buffer(glow::ARRAY_BUFFER, Some(vbo));

            let stride = (FLOATS_PER_PARTICLE * size_of::<f32>()) as i32;
            // a_pos: vec2 @ offset 0
            gl.enable_vertex_attrib_array(0);
            gl.vertex_attrib_pointer_f32(0, 2, glow::FLOAT, false, stride, 0);
            // a_radius: float @ offset 8
            gl.enable_vertex_attrib_array(1);
            gl.vertex_attrib_pointer_f32(1, 1, glow::FLOAT, false, stride, 2 * size_of::<f32>() as i32);
            // a_charge: float @ offset 12
            gl.enable_vertex_attrib_array(2);
            gl.vertex_attrib_pointer_f32(2, 1, glow::FLOAT, false, stride, 3 * size_of::<f32>() as i32);

            gl.bind_vertex_array(None);
            gl.bind_buffer(glow::ARRAY_BUFFER, None);

            let u_win_size = gl.get_uniform_location(program, "u_win_size");
            let u_coulomb = gl.get_uniform_location(program, "u_coulomb");

            gl.enable(glow::PROGRAM_POINT_SIZE);
            gl.clear_color(0.0, 0.0, 0.0, 1.0);

            Self {
                gl,
                program,
                vao,
                vbo,
                u_win_size,
                u_coulomb,
                verts: Vec::new(),
            }
        }
    }

    pub fn render(
        &mut self,
        particles: impl Iterator<Item = (f32, f32, f32, f32)>,
        coulomb_enabled: bool,
        win_size: u32,
    ) {
        self.verts.clear();
        for (x, y, r, charge) in particles {
            self.verts.extend_from_slice(&[x, y, r, charge]);
        }
        let count = (self.verts.len() / FLOATS_PER_PARTICLE) as i32;

        let gl = &self.gl;
        unsafe {
            gl.viewport(0, 0, win_size as i32, win_size as i32);
            gl.clear(glow::COLOR_BUFFER_BIT);

            if count > 0 {
                let bytes = core::slice::from_raw_parts(
                    self.verts.as_ptr() as *const u8,
                    self.verts.len() * size_of::<f32>(),
                );
                gl.bind_buffer(glow::ARRAY_BUFFER, Some(self.vbo));
                gl.buffer_data_u8_slice(glow::ARRAY_BUFFER, bytes, glow::STREAM_DRAW);

                gl.use_program(Some(self.program));
                gl.uniform_1_f32(self.u_win_size.as_ref(), win_size as f32);
                gl.uniform_1_i32(self.u_coulomb.as_ref(), coulomb_enabled as i32);

                gl.bind_vertex_array(Some(self.vao));
                gl.draw_arrays(glow::POINTS, 0, count);
                gl.bind_vertex_array(None);
            }
        }
    }
}

unsafe fn link_program(gl: &glow::Context, vs_src: &str, fs_src: &str) -> glow::Program {
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
