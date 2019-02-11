#version 430
uniform usampler2DArray sliceMap;
uniform float resolution;
uniform float sliceSize;
uniform float zOffset;


in VertexData {
    flat vec3 pos;
    vec3 normal;
} VertexIn;

layout(location = 0) out uvec4 fragColor;

void main()
{
    //discard;
    //float slice = gl_Layer*particleSize*32*4;
    uint z = uint(((VertexIn.pos.z+zOffset))/(resolution));
    //texture(sliceMap,vec3(gl_FragCoord.x,gl_FragCoord.y,gl_Layer))
    uvec4 c = uvec4(0);

    for(uint d=0;d<32;d++)
    {
        if(z!=0)
        {
            c.r |= uint(1<<d);
            z--;
        }
    }
    for(uint d=0;d<32;d++)
    {
        if(z!=0)
        {
            c.g |= uint(1<<d);
            z--;
        }
    }
    for(uint d=0;d<32;d++)
    {
        if(z!=0)
        {
            c.b |= uint(1<<d);
            z--;
        }
    }
    for(uint d=0;d<32;d++)
    {
        if(z!=0)
        {
            c.a |= uint(1<<d);
            z--;
        }
    }
    /*
    if(z<64)
    {
        c.g |= uint(1<<(z-32));
    }
    if(z<96)
    {
        c.b |= uint(1<<(z-64));
    }
    if(z<128)
    {
        c.a |= uint(1<<(z-96));
    }*/

    //fragColor = uvec4(1,2,3,4);
    fragColor = c;
}
