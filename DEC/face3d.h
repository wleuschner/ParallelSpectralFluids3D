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
            int e1;
            int e2;
            int e3;
            int e4;
        };
        int e[4];
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
