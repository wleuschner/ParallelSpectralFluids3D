#version 410

layout (location=0) in vec4 pos;

uniform mat4 pvm;

out float time;

void main()
{
    time = float(pos.w);
    gl_Position = vec4(pvm * vec4(pos.xyz,1.0));
}
