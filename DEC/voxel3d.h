#ifndef VOXEL_3D_H
#define VOXEL_3D_H
#include<glm/glm.hpp>
#include"gridenums.h"

class Voxel3D
{
public:
    Voxel3D();
    Voxel3D(unsigned int id,GridState inside);

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
};

#endif
