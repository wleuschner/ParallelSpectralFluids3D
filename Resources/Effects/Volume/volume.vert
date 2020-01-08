#version 410
layout (location=0) in vec4 pos;

void main()
{
    gl_Position = vec4(pos.x,pos.y,-1.0,1.0);
}
