#version 430
layout (location=0) in vec3 pos;
layout (location=1) in vec3 normal;

uniform float resolution;
uniform mat4 projection;

out VertexData {
    vec3 pos;
    vec3 normal;
} VertexOut;

void main()
{
    //slice = uint(position.z/(32*4*particleSize));
    VertexOut.pos = vec3(projection * vec4(pos,1.0));
    VertexOut.normal = normal;
    gl_Position = vec4(VertexOut.pos,1.0);
}
