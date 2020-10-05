#version 330
struct LightSource
{
    vec3 pos;
    vec3 ldir;
    vec3 amb;
    vec3 dif;
    vec3 spec;
};

uniform LightSource light;


uniform vec4 viewport_size;
uniform vec3 camera_position;

uniform vec3 aabb_min;
uniform vec3 aabb_max;
uniform mat4 view_mat;
uniform float step_size;
uniform float rnd;
uniform vec3 jitter;
layout(r16f) uniform sampler3D volumeTexture;

layout(location=0) out vec4 frag_colour;
in vec3 light_pos;

float rand2D(in vec2 co){
   return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453)-0.5;
}

void main()
{
    vec3 aabb_ext = vec3(aabb_max-aabb_min);
    vec3 ray_dir;
    ray_dir.xy = 2.0 * gl_FragCoord.xy / viewport_size.zw - vec2(1.0);
    ray_dir.x *= viewport_size.z/viewport_size.w;
    ray_dir.z = -1.0/tan(45.0/2.0);
    ray_dir = (vec4(ray_dir, 0) * view_mat).xyz;

    vec3 direction_inv = 1.0 / ray_dir;
    vec3 t_top = direction_inv * (aabb_max.xyz - camera_position);
    vec3 t_bottom = direction_inv * (aabb_min.xyz - camera_position);
    vec3 t_min = min(t_top, t_bottom);
    vec2 t = max(t_min.xx, t_min.yz);
    float t_0 = max(0.0, max(t.x, t.y));
    vec3 t_max = max(t_top, t_bottom);
    t = min(t_max.xx, t_max.yz);
    float t_1 = min(t.x, t.y);

    vec3 ray_start = (camera_position + ray_dir * t_0 - aabb_min.xyz) / (aabb_max.xyz - aabb_min.xyz);
    vec3 ray_stop = (camera_position + ray_dir * t_1 - aabb_min.xyz) / (aabb_max.xyz - aabb_min.xyz);

    vec3 ray = ray_stop - ray_start;
    float ray_length = length(ray);
    vec3 step_vector = step_size * ray / ray_length;

    vec3 position = ray_start;

    // Stop when the end of the volume is reached
    float transmittance = 1.0;
    vec3 color = vec3(0.0);
    float alpha = 0.0;
    int maxSteps = 1024;
    vec3 light_energy = vec3(0.0,0.0,0.0);
    float shadowthresh = -log(0.001) / 1.0;
    vec3 volVoxDims = aabb_ext/textureSize(volumeTexture,0);
    float volVoxelSize = (volVoxDims.x)*(volVoxDims.y)*(volVoxDims.z);
    while (ray_length > 0 && maxSteps>0) {
        float current_sample = (texture(volumeTexture, position).r/(volVoxelSize));
        if(current_sample>0.001)
        {
            vec3 lvec = step_size*(light_pos-position)/length(light_pos-position);
            vec3 lpos = position;
            float shadow_dist = 0;
            for(int i=0;i<128;i++)
            {
                lpos += lvec;
                float light_sample = (texture(volumeTexture, lpos).r/(volVoxelSize));
                shadow_dist+=light_sample;
                if(shadow_dist>shadowthresh) break;
            }
            float current_density = clamp(current_sample*step_size,0.0,1.0);
            vec3 shadowterm = vec3(1.0,0.5,0.0)*exp(-clamp(shadow_dist*step_size,0.0,1.0));
            vec3 absorbed_light = vec3(shadowterm*current_density);
            light_energy+=absorbed_light*transmittance;
            transmittance *= 1.0-current_density;
            if(transmittance<0.01) break;
            light_energy = clamp(light_energy,vec3(0.0),vec3(1.0));
            transmittance = clamp(transmittance,0.0,1.0);
        }

        ray_length -= step_size;
        position += step_vector;//*(1.0+0.01*rnd*vec3(0.5*cos(position.x)+0.5,0.5*sin(position.y)+0.5,0.5*cos(-position.x)+0.5));
        maxSteps--;
    }

    frag_colour = vec4(light_energy,transmittance);
}
