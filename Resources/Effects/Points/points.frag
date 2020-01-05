#version 330

layout(location=0) out vec4 frag_colour;
uniform vec4 color;
uniform float lifeTime;
in float time;

void main()
{
    frag_colour = (time/lifeTime)*color;
}
