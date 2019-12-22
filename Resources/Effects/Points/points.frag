#version 330

layout(location=0) out vec4 frag_colour;
uniform vec4 color;

void main()
{
    frag_colour = color;
}
