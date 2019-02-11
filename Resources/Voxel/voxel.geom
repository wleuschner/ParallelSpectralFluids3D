#version 430
layout(triangles) in;
layout(triangle_strip, max_vertices=3) out;
uniform float resolution;

in VertexData {
    vec3 pos;
    vec3 normal;
} VertexIn[3];

out VertexData {
    flat vec3 pos;
    vec3 normal;
} VertexOut;

uniform float zOffset;

void main()
{
    for(int i = 0;i < gl_in.length();i++)
    {
        gl_Layer = int(floor((gl_in[i].gl_Position.z+zOffset)/(32*4*resolution)));
        gl_Position = gl_in[i].gl_Position;
        VertexOut.normal = VertexIn[i].normal;
        VertexOut.pos = vec3(gl_in[i].gl_Position);
        EmitVertex();
    }
    EndPrimitive();
}

