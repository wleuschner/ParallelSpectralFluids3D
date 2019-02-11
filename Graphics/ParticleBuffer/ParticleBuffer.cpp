#include"ParticleBuffer.h"
#include<GL/glew.h>
#include<cmath>

ParticleBuffer::ParticleBuffer()
{
    glGenBuffers(1,&id);
}

ParticleBuffer::~ParticleBuffer()
{
    glDeleteBuffers(1,&id);
}

void ParticleBuffer::bind()
{
    glBindBuffer(GL_ARRAY_BUFFER,id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,id);
}

void ParticleBuffer::addParticle(Particle particle)
{
    aabb.min.x = std::fmin(aabb.min.x,particle.pos.x);
    aabb.min.y = std::fmin(aabb.min.y,particle.pos.y);
    aabb.min.z = std::fmin(aabb.min.z,particle.pos.z);
    aabb.max.x = std::fmax(aabb.max.x,particle.pos.x);
    aabb.max.y = std::fmax(aabb.max.y,particle.pos.y);
    aabb.max.z = std::fmax(aabb.max.z,particle.pos.z);
    particle.index = getNumParticles()+1;
    particles.push_back(particle);
}

void ParticleBuffer::syncGPU()
{
    bind();
    Particle* data = (Particle*)glMapBuffer(GL_ARRAY_BUFFER,GL_READ_ONLY);
    unsigned int parts = particles.size();
    particles.clear();
    for(unsigned int i=0;i<parts;i++)
    {
        particles.push_back(data[i]);
    }
    upload();
}

void ParticleBuffer::merge(const ParticleBuffer& b2)
{
    unsigned int offset = getNumParticles();
    bind();
    for(unsigned int i=0;i<b2.particles.size();i++)
    {
        Particle p = b2.particles[i];
        aabb.min = glm::min(aabb.min,glm::vec4(p.pos,0.0));
        aabb.max = glm::max(aabb.max,glm::vec4(p.pos,0.0));
        p.index = offset+i;
        addParticle(p);
    }
    upload();
}

void ParticleBuffer::upload()
{
    glBufferData(GL_ARRAY_BUFFER,particles.size()*sizeof(Particle),(void*)particles.data(),GL_DYNAMIC_DRAW);
}

void ParticleBuffer::clear()
{
    particles.clear();
}

void ParticleBuffer::updateBounds()
{
    aabb.max = glm::vec4(0.0,0.0,0.0,0.0);
    aabb.min = glm::vec4(0.0,0.0,0.0,0.0);
    for(unsigned int i=0;i<particles.size();i++)
    {
        aabb.min.x = std::fmin(aabb.min.x,particles[i].pos.x);
        aabb.min.y = std::fmin(aabb.min.y,particles[i].pos.y);
        aabb.min.z = std::fmin(aabb.min.z,particles[i].pos.z);
        aabb.max.x = std::fmax(aabb.max.x,particles[i].pos.x);
        aabb.max.y = std::fmax(aabb.max.y,particles[i].pos.y);
        aabb.max.z = std::fmax(aabb.max.z,particles[i].pos.z);
    }
}

const AABB& ParticleBuffer::getBounds()
{
    return aabb;
}

void ParticleBuffer::setBounds(const AABB& aabb)
{
    this->aabb = aabb;
}

unsigned int ParticleBuffer::getNumParticles()
{
    return particles.size();
}

std::vector<Particle>& ParticleBuffer::getParticles()
{
    return particles;
}
