#include "edge3d.h"
#include "../Geometry/mesh3d.h"
#include <cmath>

Edge3D::Edge3D()
{
    this->id = 0;
    inside = GridState::UNINITIALIZED;
}

Edge3D::Edge3D(int id,
               unsigned int v1,unsigned int v2,
               GridState inside) : Edge3D()
{
    this->inside = inside;
    this->id = id;
    this->v1 = v1;
    this->v2 = v2;
}
