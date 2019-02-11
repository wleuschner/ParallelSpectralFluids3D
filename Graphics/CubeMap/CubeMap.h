#ifndef __CUBEMAP_H_
#define __CUBEMAP_H_
#include<vector>
#include"../Texture/Texture.h"
class CubeMap
{
public:
    CubeMap();
    ~CubeMap();
    void upload_side(unsigned int side,unsigned int width,unsigned int height,void* data);
    void bind(unsigned int texUnit);
    static void unbind(unsigned int texUnit);
private:
    unsigned int id;
    unsigned int sampler;
};

#endif
