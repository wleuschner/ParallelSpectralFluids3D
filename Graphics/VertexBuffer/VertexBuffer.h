#ifndef __VERTEX_BUFFER_H
#define __VERTEX_BUFFER_H
#include<vector>
#include"../Vertex/Vertex.h"

class VertexBuffer
{
public:
    VertexBuffer();
    void bind();
    void bindBufferBase(unsigned int id);
    void upload(const std::vector<Vertex>& vertices);
private:
    unsigned int id;
};

#endif
