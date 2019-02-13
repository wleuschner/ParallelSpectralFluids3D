#include "decmesh3d.h"
#include <iostream>

DECMesh3D::DECMesh3D()
{

}

DECMesh3D::DECMesh3D(float resolution,glm::uvec3 dims,float voxelSize,glm::vec3 min)
{
    voxels.resize(dims.x*dims.y*dims.z);
    faces.resize((dims.x*dims.y)*(dims.z+1)+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))*(dims.z)); //Check if correct
    edges.resize((2*(dims.x+1)-1)*(2*(dims.y+1)-1)*(2*(dims.z+1)-1)-(dims.x+1)*(dims.y+1)*(dims.z+1));
    points.resize((dims.x+1)*(dims.y+1)*(dims.z+1));
    this->min=min;
    this->resolution = resolution;
    this->dims = dims;
    this->voxelSize = voxelSize;
    for(unsigned int z=0;z<dims.z;z++)
    {
        for(unsigned int y=0;y<dims.y;y++)
        {
            for(unsigned int x=0;x<dims.x;x++)
            {
                unsigned int v1 = (((dims.x+1)*(dims.y+1))*z)+((dims.x+1)*(y+1))+x;
                unsigned int v2 = (((dims.x+1)*(dims.y+1))*z)+((dims.x+1)*(y+1))+(x+1);
                unsigned int v3 = (((dims.x+1)*(dims.y+1))*((z+1)))+((dims.x+1)*(y+1))+(x+1);
                unsigned int v4 = (((dims.x+1)*(dims.y+1))*((z+1)))+((dims.x+1)*(y+1))+x;

                unsigned int v5 = (((dims.x+1)*(dims.y+1))*z)+((dims.x+1)*y)+x;
                unsigned int v6 = (((dims.x+1)*(dims.y+1))*z)+((dims.x+1)*y)+(x+1);
                unsigned int v7 = (((dims.x+1)*(dims.y+1))*((z+1)))+((dims.x+1)*y)+(x+1);
                unsigned int v8 = (((dims.x+1)*(dims.y+1))*((z+1)))+((dims.x+1)*y)+x;
                addVoxel(Voxel3D(z*(dims.x*dims.y)+y*dims.x+x,v1,v2,v3,v4,v5,v6,v7,v8,GridState::OUTSIDE));
            }
        }
    }
}

void DECMesh3D::addPoint(const Vertex3D& v)
{
    if(points[v.id].inside==GridState::UNINITIALIZED)
    {
        points[v.id] = v;
    }
    else if(points[v.id].inside==GridState::OUTSIDE&&
            v.inside==GridState::INSIDE)
    {
        points[v.id].inside = v.inside;
    }
}

void DECMesh3D::addEdge(const Edge3D& e)
{
    if(edges[e.id].inside==GridState::UNINITIALIZED)
    {
        edges[e.id] = e;
        addPoint(Vertex3D(e.v1,e.inside));
        addPoint(Vertex3D(e.v2,e.inside));
    }
    else if(edges[e.id].inside==GridState::OUTSIDE&&
            e.inside==GridState::INSIDE)
    {
        edges[e.id].inside = e.inside;
        addPoint(Vertex3D(e.v1,e.inside));
        addPoint(Vertex3D(e.v2,e.inside));
    }
    else
    {
        std::cout<<"EDGE CONFLICT:"<<e.id<<" "<<edges[e.id].inside<<std::endl;
    }
}

void DECMesh3D::addFace(const Face3D& f,FaceDirection direction)
{
    if(faces[f.id].inside==GridState::UNINITIALIZED)
    {
        faces[f.id] = f;
        if(direction==FaceDirection::FRONT)
        {
            faces[f.id].normal = glm::dvec3(0.0,0.0,1.0);
            unsigned int zOfs = f.id/(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            unsigned int yOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1)))/(dims.x);
            unsigned int xOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1)))%(dims.x);
            faces[f.id].e1 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs+1)*(2*dims.x+1)+xOfs;
            faces[f.id].e2 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+(xOfs+1)+dims.x;
            faces[f.id].e3 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+xOfs;
            faces[f.id].e4 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+xOfs+dims.x;
            addEdge(Edge3D(faces[f.id].e1,f.v1,f.v2,f.inside));
            addEdge(Edge3D(faces[f.id].e2,f.v2,f.v3,f.inside));
            addEdge(Edge3D(faces[f.id].e3,f.v3,f.v4,f.inside));
            addEdge(Edge3D(faces[f.id].e4,f.v4,f.v1,f.inside));

        }
        else if(direction==FaceDirection::BACK)
        {
            faces[f.id].normal = glm::dvec3(0.0,0.0,-1.0);
            unsigned int zOfs = f.id/(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            unsigned int yOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1)))/(dims.x);
            unsigned int xOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1)))%(dims.x);
            faces[f.id].e1 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(2*dims.x+1)+xOfs;
            faces[f.id].e2 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+(xOfs+1)+dims.x;
            faces[f.id].e3 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs+1)*(2*dims.x+1)+xOfs;
            faces[f.id].e4 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+(xOfs)+dims.x;
            addEdge(Edge3D(faces[f.id].e1,f.v1,f.v2,f.inside));
            addEdge(Edge3D(faces[f.id].e2,f.v2,f.v3,f.inside));
            addEdge(Edge3D(faces[f.id].e3,f.v3,f.v4,f.inside));
            addEdge(Edge3D(faces[f.id].e4,f.v4,f.v1,f.inside));

        }
        else if(direction==FaceDirection::TOP)
        {
            faces[f.id].normal = glm::dvec3(0.0,-1.0,0.0);
            unsigned int zOfs = f.id/(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            unsigned int yOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)/(2*dims.x+1);
            unsigned int xOfs = ((f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)%(2*dims.x+1));
            faces[f.id].e1 = (zOfs+1)*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(2*dims.x+1)+xOfs;
            faces[f.id].e2 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(dims.x+1)+(xOfs+1)+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            faces[f.id].e3 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(2*dims.x+1)+xOfs;
            faces[f.id].e4 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(dims.x+1)+xOfs+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            addEdge(Edge3D(faces[f.id].e1,f.v1,f.v2,f.inside));
            addEdge(Edge3D(faces[f.id].e2,f.v2,f.v3,f.inside));
            addEdge(Edge3D(faces[f.id].e3,f.v3,f.v4,f.inside));
            addEdge(Edge3D(faces[f.id].e4,f.v4,f.v1,f.inside));
        }
        else if(direction==FaceDirection::BOTTOM)
        {
            faces[f.id].normal = glm::dvec3(0.0,1.0,0.0);
            unsigned int zOfs = f.id/(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            unsigned int yOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)/(2*dims.x+1);
            unsigned int xOfs = ((f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)%(2*dims.x+1));
            faces[f.id].e1 = (zOfs)*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(2*dims.x+1)+xOfs;
            faces[f.id].e2 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(dims.x+1)+(xOfs+1)+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            faces[f.id].e3 = (zOfs+1)*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(2*dims.x+1)+xOfs;
            faces[f.id].e4 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(dims.x+1)+(xOfs)+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            addEdge(Edge3D(faces[f.id].e1,f.v1,f.v2,f.inside));
            addEdge(Edge3D(faces[f.id].e2,f.v2,f.v3,f.inside));
            addEdge(Edge3D(faces[f.id].e3,f.v3,f.v4,f.inside));
            addEdge(Edge3D(faces[f.id].e4,f.v4,f.v1,f.inside));

        }
        else if(direction==FaceDirection::LEFT)
        {
            faces[f.id].normal = glm::dvec3(1.0,0.0,0.0);
            unsigned int zOfs = f.id/(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            unsigned int yOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)/(2*dims.x+1);
            unsigned int xOfs = ((f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)%(2*dims.x+1))-dims.x;
            faces[f.id].e1 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs+1)*(dims.x+1)+xOfs+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));;
            faces[f.id].e2 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+xOfs+dims.x;
            faces[f.id].e3 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(dims.x+1)+xOfs+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            faces[f.id].e4 = (zOfs+1)*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+xOfs+dims.x;
            addEdge(Edge3D(faces[f.id].e1,f.v1,f.v2,f.inside));
            addEdge(Edge3D(faces[f.id].e2,f.v2,f.v3,f.inside));
            addEdge(Edge3D(faces[f.id].e3,f.v3,f.v4,f.inside));
            addEdge(Edge3D(faces[f.id].e4,f.v4,f.v1,f.inside));
        }
        else if(direction==FaceDirection::RIGHT)
        {
            faces[f.id].normal = glm::dvec3(-1.0,0.0,0.0);
            unsigned int zOfs = f.id/(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            unsigned int yOfs = (f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)/(2*dims.x+1);
            unsigned int xOfs = ((f.id%(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))-dims.x*dims.y)%(2*dims.x+1))-dims.x;
            faces[f.id].e1 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs+1)*(dims.x+1)+xOfs+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            faces[f.id].e2 = (zOfs+1)*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+xOfs+dims.x;
            faces[f.id].e3 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+(yOfs)*(dims.x+1)+xOfs+(2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1));
            faces[f.id].e4 = zOfs*(((2*(dims.x+1)*(dims.y+1))-(dims.x+1)-(dims.y+1))+((dims.x+1)*(dims.y+1)))+yOfs*(2*dims.x+1)+xOfs+dims.x;
            addEdge(Edge3D(faces[f.id].e1,f.v1,f.v2,f.inside));
            addEdge(Edge3D(faces[f.id].e2,f.v2,f.v3,f.inside));
            addEdge(Edge3D(faces[f.id].e3,f.v3,f.v4,f.inside));
            addEdge(Edge3D(faces[f.id].e4,f.v4,f.v1,f.inside));

        }
    }
    else if(faces[f.id].inside==GridState::OUTSIDE&&
            f.inside==GridState::INSIDE)
    {
        faces[f.id].inside = f.inside;
        addEdge(Edge3D(faces[f.id].e1,faces[f.id].v1,faces[f.id].v2,f.inside));
        addEdge(Edge3D(faces[f.id].e2,faces[f.id].v2,faces[f.id].v3,f.inside));
        addEdge(Edge3D(faces[f.id].e3,faces[f.id].v3,faces[f.id].v4,f.inside));
        addEdge(Edge3D(faces[f.id].e4,faces[f.id].v4,faces[f.id].v1,f.inside));
    }
    else {
        std::cout<<"FACE CONFLICT"<<std::endl;
    }
}

void DECMesh3D::addVoxel(const Voxel3D& v)
{
    if(voxels[v.id].inside==GridState::UNINITIALIZED)
    {
        unsigned int zOfs = v.id/(dims.x*dims.y);
        unsigned int yOfs = (v.id%(dims.x*dims.y))/dims.x;
        unsigned int xOfs = (v.id%(dims.x*dims.y))%dims.x;
        voxels[v.id] = v;
        voxels[v.id].f1 = zOfs*(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))+yOfs*(dims.x)+xOfs; //Front Face
        voxels[v.id].f2 = (zOfs+1)*(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))+yOfs*(dims.x)+xOfs; //Back Face
        voxels[v.id].f3 = zOfs*(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))+yOfs*(2*dims.x+1)+xOfs+dims.x*dims.y; //Bottom Face
        voxels[v.id].f4 = zOfs*(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))+(yOfs+1)*(2*dims.x+1)+xOfs+dims.x*dims.y; //Top Face
        voxels[v.id].f5 = zOfs*(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))+yOfs*(2*dims.x+1)+xOfs+dims.x+dims.x*dims.y; //Left Face
        voxels[v.id].f6 = zOfs*(dims.x*dims.y+2*(dims.x+1)*(dims.y+1)-(dims.x+1)-(dims.y+1))+yOfs*(2*dims.x+1)+(xOfs+1)+dims.x+dims.x*dims.y; //Right Face

        std::cout<<voxels[v.id].f1<<" "<<voxels[v.id].f2<<" "<<voxels[v.id].f3<<" "<<voxels[v.id].f4<<" "<<voxels[v.id].f5<<" "<<voxels[v.id].f6<<std::endl;
        addFace(Face3D(voxels[v.id].f1,v.v1,v.v2,v.v6,v.v5,v.inside),FaceDirection::FRONT);
        addFace(Face3D(voxels[v.id].f2,v.v8,v.v7,v.v3,v.v4,v.inside),FaceDirection::BACK);
        addFace(Face3D(voxels[v.id].f3,v.v5,v.v6,v.v7,v.v8,v.inside),FaceDirection::BOTTOM);
        addFace(Face3D(voxels[v.id].f4,v.v4,v.v3,v.v2,v.v1,v.inside),FaceDirection::TOP);
        addFace(Face3D(voxels[v.id].f5,v.v4,v.v1,v.v5,v.v8,v.inside),FaceDirection::LEFT);
        addFace(Face3D(voxels[v.id].f6,v.v2,v.v3,v.v7,v.v6,v.inside),FaceDirection::RIGHT);

    }
    else if(voxels[v.id].inside==GridState::OUTSIDE&&
            v.inside==GridState::INSIDE)
    {
        voxels[v.id].inside = v.inside;
        addFace(Face3D(voxels[v.id].f1,v.v1,v.v2,v.v6,v.v5,v.inside),FaceDirection::FRONT);
        addFace(Face3D(voxels[v.id].f2,v.v8,v.v7,v.v3,v.v4,v.inside),FaceDirection::BACK);
        addFace(Face3D(voxels[v.id].f3,v.v5,v.v6,v.v7,v.v8,v.inside),FaceDirection::BOTTOM);
        addFace(Face3D(voxels[v.id].f4,v.v4,v.v3,v.v2,v.v1,v.inside),FaceDirection::TOP);
        addFace(Face3D(voxels[v.id].f5,v.v4,v.v1,v.v5,v.v8,v.inside),FaceDirection::LEFT);
        addFace(Face3D(voxels[v.id].f6,v.v2,v.v3,v.v7,v.v6,v.inside),FaceDirection::RIGHT);
    }
    else {
        std::cout<<"CONFLICT"<<std::endl;
    }
}

PointIterator& DECMesh3D::getPointIteratorBegin()
{
    pointsBegin = points.begin();
    return pointsBegin;
}

EdgeIterator& DECMesh3D::getEdgeIteratorBegin()
{
    edgesBegin = edges.begin();
    return edgesBegin;
}

FaceIterator& DECMesh3D::getFaceIteratorBegin()
{
    facesBegin = faces.begin();
    return facesBegin;
}

VoxelIterator& DECMesh3D::getVoxelIteratorBegin()
{
    voxelsBegin = voxels.begin();
    return voxelsBegin;
}

PointIterator& DECMesh3D::getPointIteratorEnd()
{
    pointsEnd = points.end();
    return pointsEnd;
}

EdgeIterator& DECMesh3D::getEdgeIteratorEnd()
{
    edgesEnd = edges.end();
    return edgesEnd;
}

FaceIterator& DECMesh3D::getFaceIteratorEnd()
{
    facesEnd = faces.end();
    return facesEnd;
}

VoxelIterator& DECMesh3D::getVoxelIteratorEnd()
{
    voxelsEnd = voxels.end();
    return voxelsEnd;
}

glm::uvec3 DECMesh3D::getDimensions()
{
    return dims;
}

unsigned int DECMesh3D::getPointIndex(const Vertex3D& v)
{
    return v.id;
}

unsigned int DECMesh3D::getEdgeIndex(const Edge3D& e)
{
    return e.id;
}

unsigned int DECMesh3D::getFaceIndex(const Face3D& f)
{
    return f.id;
}

unsigned int DECMesh3D::getVoxelIndex(const Voxel3D& v)
{
    return v.id;
}

bool DECMesh3D::isPointInside(const glm::vec3& point)
{
    unsigned int yOfs = static_cast<unsigned int>((-min.y+point.y)/resolution)+1;
    unsigned int xOfs = static_cast<unsigned int>((-min.x+point.x)/resolution)+1;
    unsigned int zOfs = static_cast<unsigned int>((-min.x+point.z)/resolution)+1;

    return voxels[zOfs*(dims.y*dims.x)+yOfs*dims.x+xOfs].inside==GridState::INSIDE;
}

Vertex3D DECMesh3D::getPoint(unsigned int id)
{
    return points[id];
}

Edge3D DECMesh3D::getEdge(unsigned int id)
{
    return edges[id];
}

Face3D DECMesh3D::getFace(unsigned int id)
{
    return faces[id];
}

Voxel3D DECMesh3D::getVoxel(unsigned int id)
{
    return voxels[id];
}

int DECMesh3D::getPointSignum(unsigned int id,unsigned int v1)
{
    return 1;
}

int DECMesh3D::getEdgeSignum(unsigned int id,unsigned int v1,unsigned int v2)
{
    if(edges[id].v1==v1&&edges[id].v2==v2)
    {
        return 1;
    }
    return -1;
}

int DECMesh3D::getFaceSignum(unsigned int id,unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4)
{
    if((faces[id].v1==v1&&faces[id].v2==v2&&faces[id].v3==v3&&faces[id].v4==v4)||
       (faces[id].v2==v1&&faces[id].v3==v2&&faces[id].v4==v3&&faces[id].v1==v4)||
       (faces[id].v3==v1&&faces[id].v4==v2&&faces[id].v1==v3&&faces[id].v2==v4)||
       (faces[id].v4==v1&&faces[id].v1==v2&&faces[id].v2==v3&&faces[id].v3==v4))
    {
        return 1;
    }
    return -1;
}

int DECMesh3D::getVoxelSignum(unsigned int id,unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
                                   unsigned int v5,unsigned int v6,unsigned int v7,unsigned int v8)
{
    return 1;
}

unsigned int DECMesh3D::getNumPoints()
{
    return points.size();
}

unsigned int DECMesh3D::getNumEdges()
{
    return edges.size();
}

unsigned int DECMesh3D::getNumFaces()
{
    return faces.size();
}

unsigned int DECMesh3D::getNumVoxels()
{
    return voxels.size();
}
