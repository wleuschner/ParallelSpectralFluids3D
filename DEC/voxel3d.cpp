#include "voxel3d.h"

Voxel3D::Voxel3D()
{
    inside = GridState::UNINITIALIZED;
}

Voxel3D::Voxel3D(unsigned int id,
                 unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
                 unsigned int v5,unsigned int v6,unsigned int v7,unsigned int v8,
                 GridState inside) : Voxel3D()
{
    this->inside = inside;
    this->id = id;
    this->v1 = v1;
    this->v2 = v2;
    this->v3 = v3;
    this->v4 = v4;
    this->v5 = v5;
    this->v6 = v6;
    this->v7 = v7;
    this->v8 = v8;
}
