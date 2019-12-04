#ifndef FACE3D_H
#define FACE3D_H
#include<glm/glm.hpp>
#include"gridenums.h"

class Face3D
{
public:
    Face3D();
    Face3D(int id,
           unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
           GridState inside);

    int id;
    GridState inside;

    glm::dvec3 normal;

    union
    {
        struct
        {
            unsigned int e1;
            unsigned int e2;
            unsigned int e3;
            unsigned int e4;
        };
        unsigned int e[4];
    };

    union
    {
        struct
        {
            unsigned int v1;
            unsigned int v2;
            unsigned int v3;
            unsigned int v4;
        };
        unsigned int v[4];
    };
};

#endif // FACE2D_H
