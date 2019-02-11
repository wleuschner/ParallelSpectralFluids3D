#include"Light.h"

Light::Light() : Light(glm::vec3(0.0,0.0,0.0),glm::vec3(0.2,0.2,0.2),glm::vec3(0.5,0.5,0.5),glm::vec3(0.7,0.7,0.7),0.2f)
{

}

Light::Light(const glm::vec3& pos) : Light(pos,glm::vec3(0.2,0.2,0.2),glm::vec3(0.5,0.5,0.5),glm::vec3(0.7,0.7,0.7),0.2f)
{
}

Light::Light(const glm::vec3& pos,const glm::vec3& amb,const glm::vec3& diff,const glm::vec3& spec,float shininess)
{
    this->pos=pos;
    this->amb=amb;
    this->diff=diff;
    this->spec=spec;
    this->shininess=shininess;
}
