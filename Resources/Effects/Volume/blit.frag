#version 330 core

uniform vec4 viewport_size;

uniform sampler2D  texture_a;
uniform sampler2D  texture_b;
uniform sampler2D  texture_c;
uniform sampler2D  texture_d;

uniform sampler2D  texture_e;
uniform sampler2D  texture_f;
uniform sampler2D  texture_g;
uniform sampler2D  texture_h;

layout(location=0) out vec4 frag_colour;

void main()
{
    vec2 uv=(gl_FragCoord.xy-vec2(0.5))/viewport_size.zw;
    vec4 color = vec4(0.0);
    color += (0.5)*texture(texture_a,uv);
    color += (0.5/2)*texture(texture_b,uv);
    color += (0.5/4)*texture(texture_c,uv);
    color += (0.5/8)*texture(texture_d,uv);

    color += (0.5/16)*texture(texture_e,uv);
    color += (0.5/32)*texture(texture_f,uv);
    color += (0.5/64)*texture(texture_g,uv);
    color += (0.5/128)*texture(texture_h,uv);

    frag_colour = color;
}

