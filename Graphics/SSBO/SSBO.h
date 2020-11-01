#ifndef __SSBO_H
#define __SSBO_H

class SSBO
{
public:
    SSBO();
    ~SSBO();
    void clearSSBO();
    void bind(unsigned int id);
    void reserve(unsigned int num_verts);

    unsigned int id;
private:
};

#endif
