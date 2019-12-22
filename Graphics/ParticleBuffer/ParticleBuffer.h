#ifndef __PARTICLE_BUFFER_H_
#define __PARTICLE_BUFFER_H_
#include<vector>
#include"../Particle/Particle.h"
#include"../AABB/AABB.h"

class ParticleBuffer
{
public:
    ParticleBuffer();
    ~ParticleBuffer();
    void bind();
    void syncGPU();
    void reserve(unsigned int numParts);
    void clear();
    unsigned int getNumParticles();
    std::vector<Particle>& getParticles();

    unsigned int id;
private:
    AABB aabb;
    std::vector<Particle> particles;
};

#endif
