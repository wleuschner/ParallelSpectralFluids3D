#include "Particle.h"

Particle::Particle()
{
    lifeTime = 0.0;
    position = glm::dvec3(0.0);
}

Particle::Particle(double lifeTime,glm::dvec3 position)
{
    this->lifeTime = lifeTime;
    this->position = position;
}
