#include "vertex3d.h"
#include "../Geometry/mesh3d.h"

Vertex3D::Vertex3D()
{
    inside = GridState::UNINITIALIZED;
}

Vertex3D::Vertex3D(unsigned int id,GridState inside) : Vertex3D()
{
    this->inside = inside;
    this->id = id;
}
