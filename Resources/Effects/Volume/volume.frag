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
uniform usampler3D volumeTexture;

layout(location=0) out vec4 frag_colour;

void main()
{
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
    vec3 light_energy = vec3(0.0);
    while (ray_length > 0 && maxSteps>0) {
        float current_sample = texture(volumeTexture, position).r/1000.0;
        if(current_sample>0.001)
        {
            vec3 lvec = step_size*(position-light.pos);
            vec3 lpos = light.pos;
            float shadow_dist = 0;
            for(int i=0;i<128;i++)
            {
                lpos += lvec;
                float light_sample = texture(volumeTexture, lpos).r/1000.0;
                shadow_dist+=light_sample;
            }
            float current_density = current_sample*step_size;
            float shadowterm = exp(-shadow_dist*step_size);
            vec3 absorbed_light = vec3(shadowterm*current_density);
            light_energy+=absorbed_light*transmittance;
            transmittance *= 1.0-current_density;
        }

        ray_length -= step_size;
        position += step_vector;
        maxSteps--;
    }

    frag_colour = vec4(light_energy,transmittance);
}
