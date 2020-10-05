#version 410
struct LightSource
{
    vec3 pos;
    vec3 ldir;
    vec3 amb;
    vec3 dif;
    vec3 spec;
};

uniform LightSource light;
uniform mat4 view_mat;

layout (location=0) in vec4 pos;
out vec3 light_pos;

void main()
{
    light_pos = (view_mat * vec4(light.pos,1.0)).xyz;
    gl_Position = vec4(pos.x,pos.y,-1.0,1.0);
}
