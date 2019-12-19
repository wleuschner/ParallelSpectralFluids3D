#ifndef FACE3D_H
#define FACE3D_H
#include<glm/glm.hpp>
#include"gridenums.h"

class Face3D
{
public:
    Face3D();
    Face3D(int id, GridState inside);

    int id;
    GridState inside;

    glm::dvec3 center;
    glm::dvec3 normal;

    double signum;

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
            int v1;
            int v2;
        };
        int v[2];
    };
};

#endif // FACE2D_H
