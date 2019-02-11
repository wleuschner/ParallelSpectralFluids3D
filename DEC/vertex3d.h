#ifndef VERTEX3D_H
#define VERTEX3D_H
#include"gridenums.h"

class Vertex3D
{
public:
    Vertex3D();
    Vertex3D(unsigned int id,GridState inside);

    unsigned int id;
    GridState inside;
};

#endif // VERTEX2D_H
