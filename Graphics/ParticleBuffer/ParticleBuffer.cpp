#include"ParticleBuffer.h"
#include<GL/glew.h>
#include<cmath>
#include<iostream>
#include<cstdlib>

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

void ParticleBuffer::bindCompute(unsigned int idx)
{
    glBindBuffer(GL_SHADER_STORAGE_BUFFER,id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,idx,id);
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
    Particle* initData = (Particle*)malloc(num_verts*sizeof(Particle));
    memset(initData,0,num_verts*sizeof(Particle));
    glBufferData(GL_ARRAY_BUFFER,num_verts*sizeof(Particle),initData,GL_STREAM_DRAW);
    free(initData);
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
