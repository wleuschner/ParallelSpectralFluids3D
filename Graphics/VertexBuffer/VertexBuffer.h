#ifndef __VERTEX_BUFFER_H
#define __VERTEX_BUFFER_H
#include<vector>
#include"../Vertex/Vertex.h"

class VertexBuffer
{
public:
    VertexBuffer();
    ~VertexBuffer();
    void bind();
    void bindBufferBase(unsigned int id);
    void reserve(uint32_t num_verts);
    void upload(const std::vector<Vertex>& vertices);

    unsigned int id;
private:
};

#endif
