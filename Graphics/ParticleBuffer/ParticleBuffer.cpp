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
}

void ParticleBuffer::clear()
{
    glDeleteBuffers(1,&id);
    glGenBuffers(1,&id);
    particles.clear();
}

void ParticleBuffer::reserve(unsigned int num_verts)
{
    bind();
    glBufferData(GL_ARRAY_BUFFER,num_verts*sizeof(Particle),0,GL_STREAM_DRAW);
    particles.resize(num_verts);
}

void ParticleBuffer::syncGPU()
{
    bind();
    glBufferData(GL_ARRAY_BUFFER,particles.size()*sizeof(Particle),(void*)particles.data(),GL_STREAM_DRAW);
}

unsigned int ParticleBuffer::getNumParticles()
{
    return particles.size();
}

std::vector<Particle>& ParticleBuffer::getParticles()
{
    return particles;
}
