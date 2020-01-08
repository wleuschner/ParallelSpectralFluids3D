#include"Texture3D.h"
#include<GL/glew.h>
#include<cstdlib>
#include<iostream>

Texture3D::Texture3D()
{
    glGenTextures(1,&id);
    glGenSamplers(1,&sampler);
}

Texture3D::~Texture3D()
{
    destroy();
}

void Texture3D::bind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_3D,id);
    glBindSampler(texUnit,sampler);
}

void Texture3D::unbind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_3D,0);
    glBindSampler(texUnit,0);
}

void Texture3D::bindCompute(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_3D,id);
    glBindImageTexture(texUnit,id,0,GL_TRUE,0,GL_READ_WRITE,GL_R32I);
}

void Texture3D::upload(unsigned int w,unsigned int h,unsigned int d,void* data)
{
    this->width = w;
    this->height = h;
    this->depth = d;

    glTexImage3D(GL_TEXTURE_3D,0,GL_R32I,w,h,d,0,GL_RED_INTEGER,GL_INT,data);

    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);

    /*
    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glTexParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);*/
}

void Texture3D::createRenderImage(unsigned int w,unsigned int h,unsigned int d)
{

    this->width = w;
    this->height = h;
    this->depth = d;

    int result = glGetError();
    if(result!=GL_NO_ERROR)
    {
        std::cout<<"Some Error"<<result<<std::endl;
    }

    glTexImage3D(GL_TEXTURE_3D,0,GL_R32I,w,h,d,0,GL_RED_INTEGER,GL_INT,0);
    if(result!=GL_NO_ERROR)
    {
        std::cout<<"Unable to create 3D Texture "<<result<<std::endl;
    }

    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);


    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_BORDER);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_BORDER);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_BORDER);
}

void Texture3D::clearImage()
{
    glClearTexImage(id,0,GL_RED_INTEGER,GL_UNSIGNED_INT,0);
}

void Texture3D::destroy()
{
    glDeleteSamplers(1,&sampler);
    glDeleteTextures(1,&id);
}
