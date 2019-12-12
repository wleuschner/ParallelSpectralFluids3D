#include "voxel3d.h"

Voxel3D::Voxel3D()
{
    inside = GridState::UNINITIALIZED;
    this->f1 = 0;
    this->f2 = 0;
    this->f3 = 0;
    this->f4 = 0;
    this->f5 = 0;
    this->f6 = 0;
}

Voxel3D::Voxel3D(unsigned int id, GridState inside) : Voxel3D()
{
    this->inside = inside;
    this->id = id;
    this->f1 = 0;
    this->f2 = 0;
    this->f3 = 0;
    this->f4 = 0;
    this->f5 = 0;
    this->f6 = 0;
}
