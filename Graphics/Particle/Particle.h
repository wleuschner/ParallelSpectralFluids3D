#ifndef PARTICLE_H
#define PARTICLE_H
#include<glm/glm.hpp>

class Particle
{
public:
    Particle();
    Particle(double lifeTime,glm::dvec3 position);

    double lifeTime;
    glm::dvec3 position;
};

#endif
