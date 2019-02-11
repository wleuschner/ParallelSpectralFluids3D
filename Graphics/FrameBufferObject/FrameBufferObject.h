#ifndef __FRAMEBUFFEROBJECT_H
#define __FRAMEBUFFEROBJECT_H
#include"../Texture/Texture.h"
#include"../TextureArray/TextureArray.h"
#include<vector>

class FrameBufferObject
{
public:
    FrameBufferObject();
    ~FrameBufferObject();
    void bind();
    void unbind();
    void resize(unsigned int w,unsigned int h);

    void setRenderBuffer(const std::vector<unsigned int>& buffers);

    void attachColorImage(const Texture& image, unsigned int attNo);
    void attachColorArray(const TextureArray& image,unsigned int attNo);
    void attachDepthImage(const Texture& image);
    void attachDepthArray(const TextureArray& image);
    void attachStencilImage(const Texture& image);

    bool isComplete();
public:
    unsigned int id;
};

#endif
