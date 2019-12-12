#include "face3d.h"

Face3D::Face3D()
{
    this->id=0;
    inside = GridState::UNINITIALIZED;
    this->e1 = 0;
    this->e2 = 0;
    this->e3 = 0;
    this->e4 = 0;
}

Face3D::Face3D(int id, GridState inside) : Face3D()
{
    this->inside = inside;
    this->id = id;
    this->e1 = 0;
    this->e2 = 0;
    this->e3 = 0;
    this->e4 = 0;
}
