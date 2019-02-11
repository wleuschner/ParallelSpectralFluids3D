#include "TextureArray.h"
#include <GL/glew.h>

TextureArray::TextureArray()
{
    glGenTextures(1,&id);
    glGenSamplers(1,&sampler);
}

void TextureArray::createRenderArray(unsigned int width,unsigned int height,unsigned int layers)
{
    this->layers = layers;
    this->width = width;
    this->height = height;
    std::vector<unsigned int> data(width*height*layers);

    glTexStorage3D(GL_TEXTURE_2D_ARRAY,0,GL_RGBA32UI,width,height,layers);
    glTexImage3D(GL_TEXTURE_2D_ARRAY,0,GL_RGBA32UI,width,height,layers, 0, GL_RGBA_INTEGER, GL_UNSIGNED_INT, (void*)data.data());
    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
}

void TextureArray::createDepthArray(unsigned int width,unsigned int height,unsigned int layers)
{
    this->layers = layers;
    this->width = width;
    this->height = height;

    glTexStorage3D(GL_TEXTURE_2D_ARRAY,0,GL_DEPTH_COMPONENT,width,height,layers);
    glTexImage3D(GL_TEXTURE_2D_ARRAY,0,GL_DEPTH_COMPONENT,width,height,layers, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    glSamplerParameteri(sampler,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glSamplerParameteri(sampler,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
}

TextureArray::~TextureArray()
{
    glDeleteTextures(1,&id);
}

void TextureArray::bind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_2D_ARRAY,id);
    glBindSampler(texUnit,sampler);
}

void TextureArray::unbind(unsigned int texUnit)
{
    glActiveTexture(GL_TEXTURE0+texUnit);
    glBindTexture(GL_TEXTURE_2D_ARRAY,0);
    glBindSampler(texUnit,0);
}
