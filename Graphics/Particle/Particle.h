#ifndef PARTICLE_H
#define PARTICLE_H
#include<glm/glm.hpp>

class Particle
{
public:
    Particle();
    Particle(float lifeTime,glm::vec3 position);

    static void enableVertexAttribs();
    static void setVertexAttribs();

    glm::vec4 position;
};

#endif
