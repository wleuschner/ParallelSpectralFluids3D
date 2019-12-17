#ifndef EDGE3D_H
#define EDGE3D_H
#include"gridenums.h"
#include<glm/glm.hpp>

class Edge3D
{
public:
    Edge3D();
    Edge3D(int id, GridState inside);

    int id;
    GridState inside;
    glm::dvec3 center;

    unsigned int dualCount;

    union
    {
        struct
        {
            int v1;
            int v2;
        };
        int v[2];
    };

    union
    {
        struct
        {
            int f1;
            int f2;
            int f3;
            int f4;
        };
        int f[4];
    };
};

#endif // EDGE2D_H
