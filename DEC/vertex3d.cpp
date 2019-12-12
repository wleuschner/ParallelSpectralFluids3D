#include "vertex3d.h"
#include "../Geometry/mesh3d.h"

Vertex3D::Vertex3D()
{
    numDual = 0;
    inside = GridState::UNINITIALIZED;
    this->id = 0;
    this->e1 = 0;
    this->e2 = 0;
    this->e3 = 0;
    this->e4 = 0;
    this->e5 = 0;
    this->e6 = 0;
}

Vertex3D::Vertex3D(unsigned int id,GridState inside) : Vertex3D()
{
    numDual = 0;
    this->inside = inside;
    this->id = id;
    this->e1 = 0;
    this->e2 = 0;
    this->e3 = 0;
    this->e4 = 0;
    this->e5 = 0;
    this->e6 = 0;
}
