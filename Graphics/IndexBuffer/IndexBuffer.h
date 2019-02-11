#ifndef __INDEX_BUFFER_H
#define __INDEX_BUFFER_H
#include<vector>

class IndexBuffer
{
public:
    IndexBuffer();
    void bind();
    void upload(const std::vector<unsigned int>& indices);
private:
    unsigned int id;
};

#endif
