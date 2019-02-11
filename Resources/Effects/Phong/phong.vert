#version 330
struct LightSource
{
    vec3 pos;
    vec3 ldir;
    vec3 amb;
    vec3 dif;
    vec3 spec;
};

layout (location=0) in vec3 pos;
layout (location=1) in vec3 normal;
layout (location=2) in vec2 uv;

uniform mat4 pvm;
uniform mat4 modelView;
uniform mat3 normalMatrix;
uniform LightSource light;

out Data
{
    vec4 position;
    vec3 normal;
} DataOut;


void main()
{
    DataOut.position = pvm * vec4(pos,1.0);
    DataOut.normal = normalize(normalMatrix * normal);

    gl_Position = pvm * vec4(pos,1.0);
}
