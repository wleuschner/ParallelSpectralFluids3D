#ifndef __TEXTURE3D_H
#define __TEXTURE3D_H

class Texture3D
{
    friend class FrameBufferObject;
public:
    Texture3D();
    ~Texture3D();

    void bind(unsigned int texUnit);
    void bindFloatCompute(unsigned int texUnit);
    void unbind(unsigned int texUnit);
    void bindCompute(unsigned int texUnit);
    void upload(unsigned int w,unsigned int h,unsigned int d,void* data);
    void uploadFloat(unsigned int w,unsigned int h,unsigned int d,void* data);
    void createRenderImage(unsigned int w,unsigned int h,unsigned int d);
    void createFloatRenderImage(unsigned int w,unsigned int h,unsigned int d);
    void clearImage();
    void destroy();

private:
    unsigned int id;
    unsigned int sampler;
    unsigned int width;
    unsigned int height;
    unsigned int depth;
};

#endif
