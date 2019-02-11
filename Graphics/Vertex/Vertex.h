#ifndef __VERTEX_H_
#define __VERTEX_H_
#include<glm/glm.hpp>

class Vertex
{
public:
    Vertex();
    Vertex(const glm::vec3& pos);
    static void enableVertexAttribs();
    static void setVertexAttribs();

    glm::vec3 pos;
    glm::vec3 normal;
    glm::vec2 uv;
};

#endif
