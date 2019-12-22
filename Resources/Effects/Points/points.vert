#version 410
struct LightSource
{
    vec3 pos;
    vec3 ldir;
    vec3 amb;
    vec3 dif;
    vec3 spec;
};

layout (location=0) in dvec4 pos;

uniform dmat4 pvm;

void main()
{
    gl_Position = vec4(pvm * dvec4(pos.xyz,1.0));
}
