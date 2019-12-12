#include "edge3d.h"
#include "../Geometry/mesh3d.h"
#include <cmath>

Edge3D::Edge3D()
{
    this->id = 0;
    inside = GridState::UNINITIALIZED;
    this->f1 = 0;
    this->f2 = 0;
    this->f3 = 0;
    this->f4 = 0;
    this->dualCount = 0;
}

Edge3D::Edge3D(int id, GridState inside) : Edge3D()
{
    this->inside = inside;
    this->id = id;
    this->f1 = 0;
    this->f2 = 0;
    this->f3 = 0;
    this->f4 = 0;
    this->dualCount = 0;
}
