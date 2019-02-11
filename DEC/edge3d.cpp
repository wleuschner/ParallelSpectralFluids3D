#include "edge3d.h"
#include "../Geometry/mesh3d.h"
#include <cmath>

Edge3D::Edge3D()
{
    inside = GridState::UNINITIALIZED;
}

Edge3D::Edge3D(unsigned int id,
               unsigned int v1,unsigned int v2,
               GridState inside) : Edge3D()
{
    this->inside = inside;
    this->id = id;
    this->v1 = v1;
    this->v2 = v2;
}
