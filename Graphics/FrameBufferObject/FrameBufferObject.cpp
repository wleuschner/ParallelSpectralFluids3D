#include"FrameBufferObject.h"
#include<GL/glew.h>

FrameBufferObject::FrameBufferObject()
{
    glGenFramebuffers(1,&id);
}

FrameBufferObject::~FrameBufferObject()
{
    glDeleteFramebuffers(1,&id);
}

void FrameBufferObject::bind()
{
    glBindFramebuffer(GL_FRAMEBUFFER,id);
}

void FrameBufferObject::unbind()
{
    glBindFramebuffer(GL_FRAMEBUFFER,0);
}

void FrameBufferObject::resize(unsigned int w,unsigned int h)
{
    glFramebufferParameteri(GL_FRAMEBUFFER,GL_FRAMEBUFFER_DEFAULT_WIDTH,w);
    glFramebufferParameteri(GL_FRAMEBUFFER,GL_FRAMEBUFFER_DEFAULT_HEIGHT,h);
}

void FrameBufferObject::setRenderBuffer(const std::vector<unsigned int>& buffers)
{
    glDrawBuffers(buffers.size(),buffers.data());
}

void FrameBufferObject::attachColorImage(const Texture& image, unsigned int attNo)
{
    glFramebufferTexture(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0+attNo,image.id,0);
}

void FrameBufferObject::attachColorArray(const TextureArray& image,unsigned int attNo)
{
    glFramebufferTexture(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0+attNo,image.id,0);

    /*for(unsigned int i=0;i<image.layers;i++)
    {
        glFramebufferTextureLayer(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0+attNo+i,image.id,0,i);
    }*/
}

void FrameBufferObject::attachDepthImage(const Texture& image)
{
    glFramebufferTexture(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,image.id,0);
}

void FrameBufferObject::attachDepthArray(const TextureArray& image)
{
    glFramebufferTexture(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,image.id,0);
}

void FrameBufferObject::attachStencilImage(const Texture& image)
{
    glFramebufferTexture(GL_FRAMEBUFFER,GL_STENCIL_ATTACHMENT,image.id,0);
}

bool FrameBufferObject::isComplete()
{
    return glCheckFramebufferStatus(GL_FRAMEBUFFER)==GL_FRAMEBUFFER_COMPLETE;
}
