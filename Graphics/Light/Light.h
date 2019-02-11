#ifndef __LIGHT_H
#define __LIGHT_H
#include<glm/glm.hpp>

class Light
{
public:
    Light();
    Light(const glm::vec3& pos);
    Light(const glm::vec3& pos,const glm::vec3& amb,const glm::vec3& diff,const glm::vec3& spec,float shininess);

    glm::vec3 pos;
    glm::vec3 amb;
    glm::vec3 diff;
    glm::vec3 spec;
    float shininess;
};

#endif
