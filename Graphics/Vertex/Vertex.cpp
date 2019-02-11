#include"Vertex.h"
#include<GL/glew.h>

Vertex::Vertex()
{
    this->pos = glm::vec3(0.0f);
}

Vertex::Vertex(const glm::vec3& pos)
{
    this->pos = pos;
}

void Vertex::enableVertexAttribs()
{
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
}

void Vertex::setVertexAttribs()
{
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Vertex),(void*)0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(Vertex),(void*)12);
    glVertexAttribPointer(2,2,GL_FLOAT,GL_FALSE,sizeof(Vertex),(void*)24);

}
