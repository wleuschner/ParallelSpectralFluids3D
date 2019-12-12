#ifndef VERTEX3D_H
#define VERTEX3D_H
#include"gridenums.h"
#include<glm/glm.hpp>

class Vertex3D
{
public:
    Vertex3D();
    Vertex3D(unsigned int id,GridState inside);

    unsigned int id;
    glm::dvec3 center;
    GridState inside;

    unsigned int numDual;
    union
    {
        struct
        {
            int e1;
            int e2;
            int e3;
            int e4;
            int e5;
            int e6;
        };
        int e[6];
    };
};

#endif // VERTEX2D_H
