#ifndef __PARTICLE_BUFFER_H_
#define __PARTICLE_BUFFER_H_
#include<vector>
#include"../../Solver/Particle.h"
#include"../AABB/AABB.h"

class ParticleBuffer
{
public:
    ParticleBuffer();
    ~ParticleBuffer();
    void bind();
    void addParticle(Particle particle);
    void merge(const ParticleBuffer& b2);
    void syncGPU();
    void upload();
    void clear();
    void updateBounds();
    const AABB& getBounds();
    void setBounds(const AABB& aabb);
    unsigned int getNumParticles();
    std::vector<Particle>& getParticles();
private:
    unsigned int id;
    AABB aabb;
    std::vector<Particle> particles;
};

#endif
