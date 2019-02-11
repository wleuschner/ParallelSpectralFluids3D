#include"AABB.h"

AABB::AABB()
{
    max = glm::vec4(0.0,0.0,0.0,0.0);
    min = glm::vec4(0.0,0.0,0.0,0.0);
}

AABB::AABB(const glm::vec3& min,const glm::vec3& max)
{
    this->max = glm::vec4(max,0.0f);
    this->min = glm::vec4(min,0.0f);
}

const glm::vec3 AABB::getCenter() const
{
    return glm::vec3(min+0.5f*(max-min));
}

const glm::vec3 AABB::getExtent() const
{
    return glm::vec3(max-min);
}
