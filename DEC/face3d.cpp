#include "face3d.h"

Face3D::Face3D()
{
    inside = GridState::UNINITIALIZED;
}

Face3D::Face3D(int id,
               unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
               GridState inside) : Face3D()
{
    this->inside = inside;
    this->id = id;
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
    this->v4 = v4;
}
