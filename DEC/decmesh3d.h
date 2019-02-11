#ifndef DECMESH3D_H
#define DECMESH3D_H
#include<vector>
#include<glm/glm.hpp>
#include"vertex3d.h"
#include"edge3d.h"
#include"face3d.h"
#include"voxel3d.h"

typedef std::vector<Vertex3D>::iterator PointIterator;
typedef std::vector<Edge3D>::iterator EdgeIterator;
typedef std::vector<Face3D>::iterator FaceIterator;
typedef std::vector<Voxel3D>::iterator VoxelIterator;

class DECMesh3D
{
public:
    enum FaceDirection
    {
        FRONT,
        BACK,
        TOP,
        BOTTOM,
        LEFT,
        RIGHT
    };

    DECMesh3D();
    DECMesh3D(float resolution,glm::uvec3 dims,float voxelSize,glm::vec3 min);

    void addPoint(const Vertex3D& v);
    void addEdge(const Edge3D& e);
    void addFace(const Face3D& f,FaceDirection direction);
    void addVoxel(const Voxel3D& v);

    PointIterator& getPointIteratorBegin();
    EdgeIterator& getEdgeIteratorBegin();
    FaceIterator& getFaceIteratorBegin();
    VoxelIterator& getVoxelIteratorBegin();

    PointIterator& getPointIteratorEnd();
    EdgeIterator& getEdgeIteratorEnd();
    FaceIterator& getFaceIteratorEnd();
    VoxelIterator& getVoxelIteratorEnd();

    glm::uvec3 getDimensions();

    unsigned int getPointIndex(const Vertex3D& v);
    unsigned int getEdgeIndex(const Edge3D& e);
    unsigned int getFaceIndex(const Face3D& f);
    unsigned int getVoxelIndex(const Voxel3D& v);

    bool isPointInside(const glm::vec3& point);

    Vertex3D getPoint(unsigned int id);
    Edge3D getEdge(unsigned int id);
    Face3D getFace(unsigned int id);
    Voxel3D getVoxel(unsigned int id);

    int getPointSignum(unsigned int id,unsigned int v1);
    int getEdgeSignum(unsigned int id,unsigned int v1,unsigned int v2);
    int getFaceSignum(unsigned int id,unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4);
    int getVoxelSignum(unsigned int id,unsigned int v1,unsigned int v2,unsigned int v3,unsigned int v4,
                                       unsigned int v5,unsigned int v6,unsigned int v7,unsigned int v8);

    unsigned int getNumPoints();
    unsigned int getNumEdges();
    unsigned int getNumFaces();
    unsigned int getNumVoxels();
private:
    float voxelSize;
    float resolution;

    glm::vec3 min;
    glm::ivec3 dims;

    PointIterator pointsBegin;
    PointIterator pointsEnd;
    EdgeIterator edgesBegin;
    EdgeIterator edgesEnd;
    FaceIterator facesBegin;
    FaceIterator facesEnd;
    VoxelIterator voxelsBegin;
    VoxelIterator voxelsEnd;

    std::vector<Vertex3D> points;
    std::vector<Edge3D> edges;
    std::vector<Face3D> faces;
    std::vector<Voxel3D> voxels;
};

#endif
