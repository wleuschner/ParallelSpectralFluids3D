#ifndef EDGE3D_H
#define EDGE3D_H
#include"gridenums.h"

class Edge3D
{
public:
    Edge3D();
    Edge3D(unsigned int id,
           unsigned int v1,unsigned int v2,
           GridState inside);

    unsigned int id;
    GridState inside;

    union
    {
        struct
        {
            unsigned int v1;
            unsigned int v2;
        };
        unsigned int v[2];
    };
};

#endif // EDGE2D_H
