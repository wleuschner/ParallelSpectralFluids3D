#ifndef VOXEL_3D_H
#define VOXEL_3D_H
#include"gridenums.h"

class Voxel3D
{
public:
    Voxel3D();
    Voxel3D(unsigned int id,
            unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
            unsigned int v5,unsigned int v6,unsigned int v7,unsigned int v8,
            GridState inside);

    unsigned int id;
    GridState inside;

    union
    {
        struct
        {
            unsigned int f1;
            unsigned int f2;
            unsigned int f3;
            unsigned int f4;
            unsigned int f5;
            unsigned int f6;
        };
        unsigned int f[6];
    };

    union
    {
        struct
        {
            unsigned int v1;
            unsigned int v2;
            unsigned int v3;
            unsigned int v4;
            unsigned int v5;
            unsigned int v6;
            unsigned int v7;
            unsigned int v8;
        };
        unsigned int v[8];
    };
};

#endif
