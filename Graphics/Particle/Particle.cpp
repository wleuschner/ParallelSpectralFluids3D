#include "Particle.h"
#include<GL/glew.h>

Particle::Particle()
{
    position = glm::vec4(0.0);
}

Particle::Particle(float lifeTime,glm::vec3 position)
{
    this->position = glm::vec4(position,lifeTime);
}

void Particle::enableVertexAttribs()
{
    glEnableVertexAttribArray(0);
}

void Particle::setVertexAttribs()
{
    glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,sizeof(Particle),(void*)0);
}
