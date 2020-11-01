#include "SSBO.h"
#include<GL/glew.h>

SSBO::SSBO()
{
    glGenBuffers(1,&id);
}

SSBO::~SSBO()
{
    glDeleteBuffers(1,&id);
}

void SSBO::clearSSBO()
{
    int clearValue = 0;
    bind(0);
    glClearBufferData(GL_SHADER_STORAGE_BUFFER,GL_R32I,GL_RED_INTEGER,GL_INT,&clearValue);
}

void SSBO::bind(unsigned int id)
{
    glBindBuffer(GL_SHADER_STORAGE_BUFFER,id);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,id,this->id);
}

void SSBO::reserve(unsigned int num_verts)
{
    bind(0);
    glBufferData(GL_SHADER_STORAGE_BUFFER,num_verts*sizeof(unsigned int),0,GL_STREAM_DRAW);
}
