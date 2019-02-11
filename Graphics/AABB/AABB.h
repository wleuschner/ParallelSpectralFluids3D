#ifndef __AABB_H_
#define __AABB_H_
#include<glm/glm.hpp>

class AABB
{
public:
    AABB();
    AABB(const glm::vec3& min,const glm::vec3& max);
    const glm::vec3 getExtent() const;
    const glm::vec3 getCenter() const;

    glm::vec4 min;
    glm::vec4 max;
};

#endif
