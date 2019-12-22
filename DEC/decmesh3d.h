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

    enum EdgeDirection
    {
        POSZ,
        NEGZ,
        POSY,
        NEGY,
        POSX,
        NEGX,
    };

    DECMesh3D();
    DECMesh3D(float resolution,glm::uvec3 dims,float voxelSize,glm::vec3 min);

    void setPointInside(const Vertex3D& v);
    void setEdgeInside(const Edge3D& e);
    void setFaceInside(const Face3D& f);
    void setVoxelInside(const Voxel3D& v);

    int addPoint(const Vertex3D& v, unsigned x,unsigned y,unsigned z);
    int addEdge(const Edge3D& e, EdgeDirection direction, unsigned x,unsigned y,unsigned z, glm::dvec3 offset);
    int addFace(const Face3D& f,FaceDirection direction, unsigned x,unsigned y,unsigned z);
    void addVoxel(const Voxel3D& v, unsigned x,unsigned y,unsigned z);


    PointIterator& getPointIteratorBegin();
    EdgeIterator& getEdgeIteratorBegin();
    FaceIterator& getFaceIteratorBegin();
    VoxelIterator& getVoxelIteratorBegin();

    PointIterator& getPointIteratorEnd();
    EdgeIterator& getEdgeIteratorEnd();
    FaceIterator& getFaceIteratorEnd();
    VoxelIterator& getVoxelIteratorEnd();

    Voxel3D* getVoxels();
    Face3D* getFaces();
    Edge3D* getEdges();

    std::vector<unsigned int>& getSignBitString();

    glm::uvec3 getDimensions();

    unsigned int getPointIndex(const Vertex3D& v);
    unsigned int getEdgeIndex(const Edge3D& e);
    unsigned int getFaceIndex(const Face3D& f);
    unsigned int getVoxelIndex(const Voxel3D& v);

    bool isPointInside(const glm::vec3& point);

    Vertex3D getPoint(int id);
    Edge3D getEdge(int id);
    Face3D getFace(int id);
    Voxel3D getVoxel(int id);

    int getPointSignum(int id);
    int getEdgeSignum(int id);
    int getFaceSignum(int id);
    int getVoxelSignum(int id);

    int getFaceSignum(unsigned int vid,unsigned int fidx);

    unsigned int getNumPoints();
    unsigned int getNumEdges();
    unsigned int getNumFaces();
    unsigned int getNumVoxels();

    float resolution;

    unsigned int getVoxelIndex(unsigned int x,unsigned int y,unsigned z);

    unsigned int getZFaceIndex(unsigned int x,unsigned int y,unsigned z);
    unsigned int getYFaceIndex(unsigned int x,unsigned int y,unsigned z);
    unsigned int getXFaceIndex(unsigned int x,unsigned int y,unsigned z);

    unsigned int getZEdgeIndex(unsigned int x,unsigned int y,unsigned z);
    unsigned int getYEdgeIndex(unsigned int x,unsigned int y,unsigned z);
    unsigned int getXEdgeIndex(unsigned int x,unsigned int y,unsigned z);

    unsigned int getPointIndex(unsigned int x,unsigned int y,unsigned z);

    unsigned int signedIdToIndex(int id);
    int indexToSignedId(unsigned int index,int signum);

private:


    unsigned int numZFaces;
    unsigned int numYFaces;
    unsigned int numXFaces;

    unsigned int numXEdges;
    unsigned int numYEdges;
    unsigned int numZEdges;


    bool checkInternalState();

    float voxelSize;

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

    std::vector<unsigned> signBitString;
    std::vector<Vertex3D> points;
    std::vector<Edge3D> edges;
    std::vector<Face3D> faces;
    std::vector<Voxel3D> voxels;
};

#endif
