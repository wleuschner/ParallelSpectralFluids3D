#ifndef __TEXTURE_ARRAY_H_
#define __TEXTURE_ARRAY_H_
#include <vector>

class TextureArray
{
    friend class FrameBufferObject;
public:
    TextureArray();
    ~TextureArray();
    void bind(unsigned int texUnit);

    void createDepthArray(unsigned int width,unsigned int height,unsigned int layers);
    void createRenderArray(unsigned int width,unsigned int height,unsigned int layers);

    void unbind(unsigned int texUnit);
private:
    unsigned int id;
    unsigned int sampler;

    unsigned int width;
    unsigned int height;
    unsigned int layers;
    std::vector<unsigned int> layerTex;
};

#endif
