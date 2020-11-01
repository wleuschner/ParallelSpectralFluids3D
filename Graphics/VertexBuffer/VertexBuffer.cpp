#include"VertexBuffer.h"
#include<GL/glew.h>

VertexBuffer::VertexBuffer()
{
    glGenBuffers(1,&id);
}

VertexBuffer::~VertexBuffer()
{
    glDeleteBuffers(1,&id);
}

void VertexBuffer::bind()
{
    glBindBuffer(GL_ARRAY_BUFFER,id);
}

void VertexBuffer::bindBufferBase(unsigned int id)
{
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,id,this->id);
}

void VertexBuffer::reserve(uint32_t num_verts)
{
    bind();
    glBufferData(GL_ARRAY_BUFFER,num_verts*sizeof(Vertex),0,GL_STREAM_DRAW);

}

void VertexBuffer::upload(const std::vector<Vertex>& vertices)
{
    bind();
    glBufferData(GL_ARRAY_BUFFER,vertices.size()*sizeof(Vertex),(void*)vertices.data(),GL_STREAM_DRAW);
}
