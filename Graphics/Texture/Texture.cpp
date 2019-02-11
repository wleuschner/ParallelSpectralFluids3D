#include"Texture.h"
#include<GL/glew.h>

Texture::Texture()
{
    glGenTextures(1,&id);
    glGenSamplers(1,&sampler);
}

Texture::~Texture()
{
    destroy();
}

void Texture::bind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_2D,id);
    glBindSampler(texUnit,sampler);
}

void Texture::unbind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_2D,0);
    glBindSampler(texUnit,0);
}

void Texture::upload(unsigned int w,unsigned int h,void* data)
{
    this->width = w;
    this->height = h;

    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,w,h,0,GL_RGBA,GL_UNSIGNED_BYTE,data);
    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
}

void Texture::createRenderImage(unsigned int w,unsigned int h)
{
    this->width = w;
    this->height = h;

    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32UI,w,h,0,GL_RGBA_INTEGER,GL_UNSIGNED_INT,0);
    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
}

void Texture::createFloatRenderImage(unsigned int w,unsigned int h)
{
    this->width = w;
    this->height = h;

    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F,w,h,0,GL_RGBA,GL_FLOAT,0);
    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
}

void Texture::createDepthImage(unsigned int w,unsigned int h)
{
    this->width = w;
    this->height = h;

    glTexImage2D(GL_TEXTURE_2D,0,GL_DEPTH_COMPONENT32F,w,h,0,GL_DEPTH_COMPONENT,GL_FLOAT,0);
    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
}

void Texture::destroy()
{
    glDeleteSamplers(1,&sampler);
    glDeleteTextures(1,&id);
}
