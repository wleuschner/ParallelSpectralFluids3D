#version 330
struct LightSource
{
    vec3 pos;
    vec3 ldir;
    vec3 amb;
    vec3 dif;
    vec3 spec;
};

in Data
{
    vec4 position;
    vec3 normal;
} DataIn;

uniform LightSource light;

layout(location=0) out vec4 frag_colour;

void main()
{
    vec4 diff;
    diff = vec4(0.2,0.2,0.2,1.0);

    vec4 spec = vec4(0.0);
    float intensity = max(dot(DataIn.normal,light.ldir),0.0);
    if(intensity>0.0)
    {
        vec3 p = vec3(normalize(DataIn.position));
        vec3 e = normalize(-p);
        vec3 h = normalize(light.ldir + e);
        float intSpec = max(dot(h,DataIn.normal),0.0);
        spec = vec4(0.5,0.5,0.5,1.0) * pow(intSpec,1.2);
    }
    //frag_colour = vec4(1.0,0.0,0.0,1.0);
    frag_colour = vec4(max(intensity*diff/*+spec*/,diff).xyz,0.5);
}
