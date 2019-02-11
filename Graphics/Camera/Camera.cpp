#include"Camera.h"
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtx/rotate_vector.hpp>

Camera::Camera(glm::vec3 pos)
{
    this->eye = pos;
    this->up = glm::vec3(0.0,1.0,0.0);
    this->center = glm::vec3(0.0,0.0,1.0);
}

void Camera::translate(const glm::vec3& d)
{
    eye += d;
}

void Camera::rotate(float a,const glm::vec3& axis)
{
    center = glm::rotate(center,a,axis);
}

glm::mat4 Camera::getView()
{
    return glm::lookAt(eye,eye+center,up);
}

glm::mat4 Camera::getRotMat()
{
    return glm::mat4(glm::mat3(getView()));
}

glm::vec3 Camera::getPosition()
{
    return eye;
}

glm::vec3 Camera::getForwardVec()
{
    return center;
}

glm::vec3 Camera::getUpVector()
{
    return up;
}

glm::vec3 Camera::getStrafeVec()
{
    return glm::cross(center,up);
}
