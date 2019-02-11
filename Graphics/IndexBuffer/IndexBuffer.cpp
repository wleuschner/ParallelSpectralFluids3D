#include"IndexBuffer.h"
#include<GL/glew.h>

IndexBuffer::IndexBuffer()
{
    glGenBuffers(1,&id);
}

void IndexBuffer::bind()
{
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,id);
}

void IndexBuffer::upload(const std::vector<unsigned int>& indices)
{
    bind();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,indices.size()*sizeof(unsigned int),(void*)indices.data(),GL_STATIC_DRAW);
}
