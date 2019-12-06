#ifndef VOXEL_3D_H
#define VOXEL_3D_H
#include<glm/glm.hpp>
#include"gridenums.h"

class Voxel3D
{
public:
    Voxel3D();
    Voxel3D(unsigned int id,
            unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
            unsigned int v5,unsigned int v6,unsigned int v7,unsigned int v8,
            GridState inside);

    int id;
    GridState inside;
    glm::dvec3 center;

    union
    {
        struct
        {
            int f1;
            int f2;
            int f3;
            int f4;
            int f5;
            int f6;
        };
        int f[6];
    };

    union
    {
        struct
        {
            int v1;
            int v2;
            int v3;
            int v4;
            int v5;
            int v6;
            int v7;
            int v8;
        };
        int v[8];
    };
};

#endif
