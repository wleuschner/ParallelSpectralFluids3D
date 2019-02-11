#include"CubeMap.h"
#include<GL/glew.h>

CubeMap::CubeMap()
{
    glGenTextures(1,&id);
    glGenSamplers(1,&sampler);
    glSamplerParameteri(sampler, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glSamplerParameteri(sampler, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glSamplerParameteri(sampler, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}

CubeMap::~CubeMap()
{
    glDeleteTextures(1,&id);
    glDeleteSamplers(1,&sampler);
}

void CubeMap::bind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_CUBE_MAP,id);
    glBindSampler(texUnit,sampler);
}

void CubeMap::unbind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_CUBE_MAP,0);
    glBindSampler(texUnit,0);
}

void CubeMap::upload_side(unsigned int side,unsigned int width,unsigned int height,void* data)
{
    glTexImage2D(side,0,GL_RGBA,width,height,0,GL_BGRA,GL_UNSIGNED_BYTE,data);
}
