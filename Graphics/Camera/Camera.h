#ifndef __CAMERA_H_
#define __CAMERA_H_
#include<glm/glm.hpp>

class Camera
{
public:
    Camera(glm::vec3 pos=glm::vec3(0.0,0.0,0.0));

    void translate(const glm::vec3& d);
    void rotate(float a,const glm::vec3& axis);
    glm::mat4 getView();
    glm::mat4 getRotMat();
    glm::vec3 getPosition();
    glm::vec3 getForwardVec();
    glm::vec3 getUpVector();
    glm::vec3 getStrafeVec();
private:
    glm::vec3 up;
    glm::vec3 center;
    glm::vec3 eye;
};

#endif
