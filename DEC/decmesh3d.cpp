#include "decmesh3d.h"
#include <iostream>
#include <cstdlib>

#include <algorithm>

DECMesh3D::DECMesh3D()
{

}

unsigned int DECMesh3D::getVoxelIndex(unsigned int x,unsigned int y,unsigned z)
{
    return z*(dims.x*dims.y)+y*dims.x+x;
}

unsigned int DECMesh3D::getZFaceIndex(unsigned int x,unsigned int y,unsigned z)
{
    return (z*(dims.x*dims.y)+y*dims.x+x);
}

unsigned int DECMesh3D::getYFaceIndex(unsigned int x,unsigned int y,unsigned z)
{
    return (z*((dims.y+1)*(dims.x))+y*(dims.x)+x)+numZFaces;
}

unsigned int DECMesh3D::getXFaceIndex(unsigned int x,unsigned int y,unsigned z)
{
    return (z*((dims.y)*(dims.x+1))+y*(dims.x+1)+x+numZFaces+numYFaces);
}

unsigned int DECMesh3D::getZEdgeIndex(unsigned int x,unsigned int y,unsigned z)
{
    return z*(dims.x+1)*(dims.y+1)+y*(dims.x+1)+x;
}

unsigned int DECMesh3D::getYEdgeIndex(unsigned int x,unsigned int y,unsigned z)
{
    return z*(dims.x+1)*(dims.y)+y*(dims.x+1)+x+numZEdges;
}

unsigned int DECMesh3D::getXEdgeIndex(unsigned int x,unsigned int y,unsigned z)
{
    return z*(dims.x)*(dims.y+1)+y*(dims.x)+x+numZEdges+numYEdges;
}

unsigned int DECMesh3D::getPointIndex(unsigned int x,unsigned int y,unsigned z)
{
    return z*((dims.y+1)*(dims.x+1))+y*(dims.x+1)+x;
}

unsigned int DECMesh3D::signedIdToIndex(int id)
{
    return labs(id)-1;
}

int DECMesh3D::indexToSignedId(unsigned int index,int signum)
{
    return signum*(index+1);
}

DECMesh3D::DECMesh3D(float resolution,glm::uvec3 dims,float voxelSize,glm::vec3 min)
{
    numXFaces = (dims.x+1)*(dims.y*dims.z);
    numYFaces = (dims.y+1)*(dims.x*dims.z);
    numZFaces = (dims.z+1)*(dims.x*dims.y);

    numXEdges = (dims.x)*(dims.y+1)*(dims.z+1);
    numYEdges = (dims.x+1)*(dims.y)*(dims.z+1);
    numZEdges = (dims.x+1)*(dims.y+1)*(dims.z);

    voxels.resize(dims.x*dims.y*dims.z);
    faces.resize(numXFaces+numYFaces+numZFaces); //Check if correct
    //edges.resize((2*(dims.x+1)-1)*(2*(dims.y+1)-1)*(2*(dims.z+1)-1)-(dims.x+1)*(dims.y+1)*(dims.z+1));
    edges.resize(numXEdges+numYEdges+numZEdges);
    points.resize((dims.x+1)*(dims.y+1)*(dims.z+1));
    this->min=min;
    this->resolution = resolution;
    this->dims = dims;
    this->voxelSize = voxelSize;


    //Add Voxel
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z;z++)
    {
        for(unsigned int y=0;y<dims.y;y++)
        {
            for(unsigned int x=0;x<dims.x;x++)
            {
                unsigned int vid = z*(dims.x*dims.y)+y*dims.x+x;
                assert(vid<getNumVoxels());

                unsigned int zOfs = z;
                unsigned int yOfs = y;
                unsigned int xOfs = x;

                voxels[vid].id = vid;
                voxels[vid].inside = GridState::OUTSIDE;
                voxels[vid].center = glm::dvec3(min.x+x*resolution+resolution/2,min.y+y*resolution+resolution/2,min.z+z*resolution+resolution/2);

                voxels[vid].f1 = indexToSignedId(getZFaceIndex(x,y,z),1);//(z*((dims.z+1)*(dims.x*dims.y)+(dims.y+1)*(dims.x*dims.z)+(dims.x+1)*(dims.y*dims.z))+y*(dims.x)+x+1); //Front Face
                voxels[vid].f2 = indexToSignedId(getZFaceIndex(x,y,z+1),-1); //Back Face
                voxels[vid].f3 = indexToSignedId(getYFaceIndex(x,y,z),1); //Bottom Face
                voxels[vid].f4 = indexToSignedId(getYFaceIndex(x,y+1,z),-1); //Top Face
                voxels[vid].f5 = indexToSignedId(getXFaceIndex(x,y,z),1); //Left Face
                voxels[vid].f6 = indexToSignedId(getXFaceIndex(x+1,y,z),-1); //Right Face

                //Dual
                faces[getZFaceIndex(x,y,z)].v1=vid;
                faces[getZFaceIndex(x,y,z+1)].v2=vid;
                faces[getYFaceIndex(x,y,z)].v1=vid;
                faces[getYFaceIndex(x,y+1,z)].v2=vid;
                faces[getXFaceIndex(x,y,z)].v1=vid;
                faces[getXFaceIndex(x+1,y,z)].v2=vid;

                assert(labs(voxels[vid].f1)!=labs(voxels[vid].f2) &&
                       labs(voxels[vid].f1)!=labs(voxels[vid].f3) &&
                       labs(voxels[vid].f1)!=labs(voxels[vid].f4) &&
                       labs(voxels[vid].f1)!=labs(voxels[vid].f5) &&
                       labs(voxels[vid].f1)!=labs(voxels[vid].f6) &&
                       labs(voxels[vid].f2)!=labs(voxels[vid].f3) &&
                       labs(voxels[vid].f2)!=labs(voxels[vid].f4) &&
                       labs(voxels[vid].f2)!=labs(voxels[vid].f5) &&
                       labs(voxels[vid].f2)!=labs(voxels[vid].f6) &&
                       labs(voxels[vid].f3)!=labs(voxels[vid].f4) &&
                       labs(voxels[vid].f3)!=labs(voxels[vid].f5) &&
                       labs(voxels[vid].f3)!=labs(voxels[vid].f6) &&
                       labs(voxels[vid].f4)!=labs(voxels[vid].f5) &&
                       labs(voxels[vid].f4)!=labs(voxels[vid].f6) &&
                       labs(voxels[vid].f5)!=labs(voxels[vid].f6));

                assert(labs(voxels[vid].f1)-1<getNumFaces());
                assert(labs(voxels[vid].f2)-1<getNumFaces());
                assert(labs(voxels[vid].f3)-1<getNumFaces());
                assert(labs(voxels[vid].f4)-1<getNumFaces());
                assert(labs(voxels[vid].f5)-1<getNumFaces());
                assert(labs(voxels[vid].f6)-1<getNumFaces());
            }
        }
    }

    //Add Front Faces
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z+1;z++)
    {
        for(unsigned int y=0;y<dims.y;y++)
        {
            for(unsigned int x=0;x<dims.x;x++)
            {
                unsigned int faceFrontIdx = getZFaceIndex(x,y,z);
                assert(faceFrontIdx<getNumFaces());
                assert(faces[faceFrontIdx].id==0);
                unsigned int zOfs = z;
                unsigned int yOfs = y;
                unsigned int xOfs = x;
                //Front Face Prolog
                assert(faces[faceFrontIdx].id==0);

                faces[faceFrontIdx].normal = glm::dvec3(0.0,0.0,1.0);
                faces[faceFrontIdx].center = glm::dvec3(min.x+x*resolution+resolution/2,min.y+y*resolution+resolution/2,min.z+z*resolution);

                faces[faceFrontIdx].id = faceFrontIdx;
                faces[faceFrontIdx].inside = GridState::OUTSIDE;

                assert(faces[faceFrontIdx].e1 == 0 || faces[faceFrontIdx].e1 == indexToSignedId(getXEdgeIndex(x,y+1,z),-1));
                assert(faces[faceFrontIdx].e2 == 0 || faces[faceFrontIdx].e2 == indexToSignedId(getYEdgeIndex(x+1,y,z),-1));
                assert(faces[faceFrontIdx].e3 == 0 || faces[faceFrontIdx].e3 == indexToSignedId(getXEdgeIndex(x,y,z),1));
                assert(faces[faceFrontIdx].e4 == 0 || faces[faceFrontIdx].e4 == indexToSignedId(getYEdgeIndex(x,y,z),1));

                faces[faceFrontIdx].e1 = indexToSignedId(getXEdgeIndex(x,y+1,z),-1);
                faces[faceFrontIdx].e2 = indexToSignedId(getYEdgeIndex(x+1,y,z),-1);
                faces[faceFrontIdx].e3 = indexToSignedId(getXEdgeIndex(x,y,z),1);
                faces[faceFrontIdx].e4 = indexToSignedId(getYEdgeIndex(x,y,z),1);

                //Dual
                /*
                edges[getXEdgeIndex(x,y+1,z)].f[edges[getXEdgeIndex(x,y+1,z)].dualCount] = indexToSignedId(faceFrontIdx,1);
                edges[getYEdgeIndex(x+1,y,z)].f[edges[getYEdgeIndex(x+1,y,z)].dualCount] = indexToSignedId(faceFrontIdx,1);
                edges[getXEdgeIndex(x,y,z)].f[edges[getXEdgeIndex(x,y,z)].dualCount] = indexToSignedId(faceFrontIdx,1);
                edges[getYEdgeIndex(x,y,z)].f[edges[getYEdgeIndex(x,y,z)].dualCount] = indexToSignedId(faceFrontIdx,1);
*/
                edges[getXEdgeIndex(x,y+1,z)].f[0] = indexToSignedId(faceFrontIdx,-1);
                edges[getYEdgeIndex(x+1,y,z)].f[2] = indexToSignedId(faceFrontIdx,-1);
                edges[getXEdgeIndex(x,y,z)].f[1] = indexToSignedId(faceFrontIdx,1);
                edges[getYEdgeIndex(x,y,z)].f[3] = indexToSignedId(faceFrontIdx,1);

                /*edges[getXEdgeIndex(x,y+1,z)].dualCount++;
                edges[getYEdgeIndex(x+1,y,z)].dualCount++;
                edges[getXEdgeIndex(x,y,z)].dualCount++;
                edges[getYEdgeIndex(x,y,z)].dualCount++;*/

                /*
                edges[getXEdgeIndex(x,y+1,z)].f1 = indexToSignedId(faceFrontIdx,-1);
                edges[getYEdgeIndex(x+1,y,z)].f2 = indexToSignedId(faceFrontIdx,-1);
                edges[getXEdgeIndex(x,y,z)].f3 = indexToSignedId(faceFrontIdx,1);
                edges[getYEdgeIndex(x,y,z)].f4 = indexToSignedId(faceFrontIdx,1);*/

                assert(labs(faces[faceFrontIdx].e1)!=labs(faces[faceFrontIdx].e2) &&
                       labs(faces[faceFrontIdx].e1)!=labs(faces[faceFrontIdx].e3) &&
                       labs(faces[faceFrontIdx].e1)!=labs(faces[faceFrontIdx].e4) &&
                       labs(faces[faceFrontIdx].e2)!=labs(faces[faceFrontIdx].e3) &&
                       labs(faces[faceFrontIdx].e2)!=labs(faces[faceFrontIdx].e4) &&
                       labs(faces[faceFrontIdx].e3)!=labs(faces[faceFrontIdx].e4));
                assert(labs(faces[faceFrontIdx].e1)-1<getNumEdges());
                assert(labs(faces[faceFrontIdx].e2)-1<getNumEdges());
                assert(labs(faces[faceFrontIdx].e3)-1<getNumEdges());
                assert(labs(faces[faceFrontIdx].e4)-1<getNumEdges());

            }
        }
    }

    //Add Bottom Faces
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z;z++)
    {
        for(unsigned int y=0;y<dims.y+1;y++)
        {
            for(unsigned int x=0;x<dims.x;x++)
            {
                unsigned int faceBottomIdx = getYFaceIndex(x,y,z);
                assert(faceBottomIdx<getNumFaces());
                assert(faces[faceBottomIdx].id==0);
                unsigned int zOfs = z;
                unsigned int yOfs = y;
                unsigned int xOfs = x;

                //Bottom Face Prolog
                assert(faces[faceBottomIdx].id==0);

                faces[faceBottomIdx].normal = glm::dvec3(0.0,1.0,0.0);
                faces[faceBottomIdx].center = glm::dvec3(min.x+x*resolution+resolution/2,min.y+y*resolution,min.z+z*resolution+resolution/2);


                faces[faceBottomIdx].id = faceBottomIdx;
                faces[faceBottomIdx].inside = GridState::OUTSIDE;

                assert(faces[faceBottomIdx].e3 == 0 || faces[faceBottomIdx].e3 == indexToSignedId(getXEdgeIndex(x,y,z+1),1));
                assert(faces[faceBottomIdx].e2 == 0 || faces[faceBottomIdx].e2 == indexToSignedId(getZEdgeIndex(x+1,y,z),1));
                assert(faces[faceBottomIdx].e1 == 0 || faces[faceBottomIdx].e1 == indexToSignedId(getXEdgeIndex(x,y,z),-1));
                assert(faces[faceBottomIdx].e4 == 0 || faces[faceBottomIdx].e4 == indexToSignedId(getZEdgeIndex(x,y,z),-1));

                faces[faceBottomIdx].e3 = indexToSignedId(getXEdgeIndex(x,y,z+1),1);
                faces[faceBottomIdx].e2 = indexToSignedId(getZEdgeIndex(x+1,y,z),1);
                faces[faceBottomIdx].e1 = indexToSignedId(getXEdgeIndex(x,y,z),-1);
                faces[faceBottomIdx].e4 = indexToSignedId(getZEdgeIndex(x,y,z),-1);

                //Dual
/*
                edges[getXEdgeIndex(x,y,z+1)].f[edges[getXEdgeIndex(x,y,z+1)].dualCount] = indexToSignedId(faceBottomIdx,1);
                edges[getZEdgeIndex(x+1,y,z)].f[edges[getZEdgeIndex(x+1,y,z)].dualCount] = indexToSignedId(faceBottomIdx,1);
                edges[getXEdgeIndex(x,y,z)].f[edges[getXEdgeIndex(x,y,z)].dualCount] = indexToSignedId(faceBottomIdx,1);
                edges[getZEdgeIndex(x,y,z)].f[edges[getZEdgeIndex(x,y,z)].dualCount] = indexToSignedId(faceBottomIdx,1);
*/
                edges[getXEdgeIndex(x,y,z+1)].f[2] = indexToSignedId(faceBottomIdx,1);
                edges[getZEdgeIndex(x+1,y,z)].f[0] = indexToSignedId(faceBottomIdx,1);
                edges[getXEdgeIndex(x,y,z)].f[3] = indexToSignedId(faceBottomIdx,-1);
                edges[getZEdgeIndex(x,y,z)].f[1] = indexToSignedId(faceBottomIdx,-1);

                /*
                edges[getXEdgeIndex(x,y,z+1)].dualCount++;
                edges[getZEdgeIndex(x+1,y,z)].dualCount++;
                edges[getXEdgeIndex(x,y,z)].dualCount++;
                edges[getZEdgeIndex(x,y,z)].dualCount++;*/

                /*
                edges[getXEdgeIndex(x,y,z+1)].f3 = indexToSignedId(faceBottomIdx,1);
                edges[getZEdgeIndex(x+1,y,z)].f2 = indexToSignedId(faceBottomIdx,1);
                edges[getXEdgeIndex(x,y,z)].f1 = indexToSignedId(faceBottomIdx,-1);
                edges[getZEdgeIndex(x,y,z)].f4 = indexToSignedId(faceBottomIdx,-1);*/

                assert(labs(faces[faceBottomIdx].e1)!=labs(faces[faceBottomIdx].e2) &&
                       labs(faces[faceBottomIdx].e1)!=labs(faces[faceBottomIdx].e3) &&
                       labs(faces[faceBottomIdx].e1)!=labs(faces[faceBottomIdx].e4) &&
                       labs(faces[faceBottomIdx].e2)!=labs(faces[faceBottomIdx].e3) &&
                       labs(faces[faceBottomIdx].e2)!=labs(faces[faceBottomIdx].e4) &&
                       labs(faces[faceBottomIdx].e3)!=labs(faces[faceBottomIdx].e4));
                assert(labs(faces[faceBottomIdx].e1)-1<getNumEdges());
                assert(labs(faces[faceBottomIdx].e2)-1<getNumEdges());
                assert(labs(faces[faceBottomIdx].e3)-1<getNumEdges());
                assert(labs(faces[faceBottomIdx].e4)-1<getNumEdges());
            }
        }
    }

    //Add Left Faces
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z;z++)
    {
        for(unsigned int y=0;y<dims.y;y++)
        {
            for(unsigned int x=0;x<dims.x+1;x++)
            {
                unsigned int faceLeftIdx = getXFaceIndex(x,y,z);
                assert(faceLeftIdx<getNumFaces());
                assert(faces[faceLeftIdx].id==0);
                unsigned int zOfs = z;
                unsigned int yOfs = y;
                unsigned int xOfs = x;

                //Left Face Prolog
                faces[faceLeftIdx].normal = glm::dvec3(1.0,0.0,0.0);
                faces[faceLeftIdx].center = glm::dvec3(min.x+x*resolution,min.y+y*resolution+resolution/2,min.z+z*resolution+resolution/2);


                faces[faceLeftIdx].id = faceLeftIdx;
                faces[faceLeftIdx].inside = GridState::OUTSIDE;

                assert(faces[faceLeftIdx].e1 == 0 || faces[faceLeftIdx].e1 == indexToSignedId(getZEdgeIndex(x,y+1,z),1));
                assert(faces[faceLeftIdx].e2 == 0 || faces[faceLeftIdx].e2 == indexToSignedId(getYEdgeIndex(x,y,z),-1));
                assert(faces[faceLeftIdx].e3 == 0 || faces[faceLeftIdx].e3 == indexToSignedId(getZEdgeIndex(x,y,z),-1));
                assert(faces[faceLeftIdx].e4 == 0 || faces[faceLeftIdx].e4 == indexToSignedId(getYEdgeIndex(x,y,z+1),1));

                faces[faceLeftIdx].e1 = indexToSignedId(getZEdgeIndex(x,y+1,z),1);
                faces[faceLeftIdx].e2 = indexToSignedId(getYEdgeIndex(x,y,z),-1);
                faces[faceLeftIdx].e3 = indexToSignedId(getZEdgeIndex(x,y,z),-1);
                faces[faceLeftIdx].e4 = indexToSignedId(getYEdgeIndex(x,y,z+1),1);

                //Dual

                edges[getZEdgeIndex(x,y+1,z)].f[2] = indexToSignedId(faceLeftIdx,1);
                edges[getYEdgeIndex(x,y,z)].f[0] = indexToSignedId(faceLeftIdx,-1);
                edges[getZEdgeIndex(x,y,z)].f[3] = indexToSignedId(faceLeftIdx,-1);
                edges[getYEdgeIndex(x,y,z+1)].f[1] = indexToSignedId(faceLeftIdx,1);

                /*edges[getZEdgeIndex(x,y+1,z)].dualCount++;
                edges[getYEdgeIndex(x,y,z)].dualCount++;
                edges[getZEdgeIndex(x,y,z)].dualCount++;
                edges[getYEdgeIndex(x,y,z+1)].dualCount++;*/

                /*
                edges[getZEdgeIndex(x,y+1,z)].f1 = indexToSignedId(faceLeftIdx,1);
                edges[getYEdgeIndex(x,y,z)].f2 = indexToSignedId(faceLeftIdx,-1);
                edges[getZEdgeIndex(x,y,z)].f3 = indexToSignedId(faceLeftIdx,-1);
                edges[getYEdgeIndex(x,y,z+1)].f4 = indexToSignedId(faceLeftIdx,1);*/

                assert(labs(faces[faceLeftIdx].e1)!=labs(faces[faceLeftIdx].e2) &&
                       labs(faces[faceLeftIdx].e1)!=labs(faces[faceLeftIdx].e3) &&
                       labs(faces[faceLeftIdx].e1)!=labs(faces[faceLeftIdx].e4) &&
                       labs(faces[faceLeftIdx].e2)!=labs(faces[faceLeftIdx].e3) &&
                       labs(faces[faceLeftIdx].e2)!=labs(faces[faceLeftIdx].e4) &&
                       labs(faces[faceLeftIdx].e3)!=labs(faces[faceLeftIdx].e4));


                assert(labs(faces[faceLeftIdx].e1)-1<getNumEdges());
                assert(labs(faces[faceLeftIdx].e2)-1<getNumEdges());
                assert(labs(faces[faceLeftIdx].e3)-1<getNumEdges());
                assert(labs(faces[faceLeftIdx].e4)-1<getNumEdges());
            }
        }
    }

    //Add Bottom Left Edges
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z;z++)
    {
        for(unsigned int y=0;y<dims.y+1;y++)
        {
            for(unsigned int x=0;x<dims.x+1;x++)
            {
                unsigned int edgeBLIdx = getZEdgeIndex(x,y,z);
                assert(edgeBLIdx<getNumEdges());
                assert(edges[edgeBLIdx].id==0);

                edges[edgeBLIdx].center = glm::dvec3(min.x+x*resolution,min.y+y*resolution,min.z+z*resolution+resolution/2);
                edges[edgeBLIdx].inside = GridState::OUTSIDE;

                edges[edgeBLIdx].id = edgeBLIdx;
                edges[edgeBLIdx].v2 = getPointIndex(x,y,z);
                edges[edgeBLIdx].v1 = getPointIndex(x,y,z+1);

                points[getPointIndex(x,y,z)].e2 = indexToSignedId(edgeBLIdx,1);
                points[getPointIndex(x,y,z+1)].e1 = indexToSignedId(edgeBLIdx,-1);

                assert(edges[edgeBLIdx].v1<getNumPoints());
                assert(edges[edgeBLIdx].v2<getNumPoints());
            }
        }
    }

    //Add Bottom Front Edges
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z+1;z++)
    {
        for(unsigned int y=0;y<dims.y+1;y++)
        {
            for(unsigned int x=0;x<dims.x;x++)
            {
                unsigned int edgeBFIdx = getXEdgeIndex(x,y,z);
                assert(edgeBFIdx<getNumEdges());
                assert(edges[edgeBFIdx].id==0);

                edges[edgeBFIdx].center = glm::dvec3(min.x+x*resolution+resolution/2,min.y+y*resolution,min.z+z*resolution);
                edges[edgeBFIdx].inside = GridState::OUTSIDE;

                edges[edgeBFIdx].id = edgeBFIdx;
                edges[edgeBFIdx].v2 = getPointIndex(x,y,z);
                edges[edgeBFIdx].v1 = getPointIndex(x+1,y,z);

                points[getPointIndex(x,y,z)].e6 = indexToSignedId(edgeBFIdx,1);
                points[getPointIndex(x+1,y,z)].e5 = indexToSignedId(edgeBFIdx,-1);

                assert(edges[edgeBFIdx].v1<getNumPoints());
                assert(edges[edgeBFIdx].v2<getNumPoints());
            }
        }
    }

    //Add Left Front Edges
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z+1;z++)
    {
        for(unsigned int y=0;y<dims.y;y++)
        {
            for(unsigned int x=0;x<dims.x+1;x++)
            {
                unsigned int edgeLFIdx = getYEdgeIndex(x,y,z);
                assert(edgeLFIdx<getNumEdges());
                assert(edges[edgeLFIdx].id==0);

                edges[edgeLFIdx].center = glm::dvec3(min.x+x*resolution,min.y+y*resolution+resolution/2,min.z+z*resolution);
                edges[edgeLFIdx].inside = GridState::OUTSIDE;

                edges[edgeLFIdx].id = edgeLFIdx;
                edges[edgeLFIdx].v2 = getPointIndex(x,y,z);
                edges[edgeLFIdx].v1 = getPointIndex(x,y+1,z);

                points[getPointIndex(x,y,z)].e4 = indexToSignedId(edgeLFIdx,1);
                points[getPointIndex(x,y+1,z)].e3 = indexToSignedId(edgeLFIdx,-1);

                assert(edges[edgeLFIdx].v1<getNumPoints());
                assert(edges[edgeLFIdx].v2<getNumPoints());
            }
        }
    }

    //Add Points
    //#pragma omp parallel for
    for(unsigned int z=0;z<dims.z+1;z++)
    {
        for(unsigned int y=0;y<dims.y+1;y++)
        {
            for(unsigned int x=0;x<dims.x+1;x++)
            {
                unsigned int vid = getPointIndex(x,y,z);
                assert(vid<getNumPoints());

                points[vid].center = glm::dvec3(min.x+x*resolution,min.y+y*resolution,min.z+z*resolution);
                points[vid].id = vid;
                points[vid].inside = GridState::INSIDE;

                /*assert(labs(points[vid].e1)!=labs(points[vid].e2) &&
                       labs(points[vid].e1)!=labs(points[vid].e3) &&
                       labs(points[vid].e1)!=labs(points[vid].e4) &&
                       labs(points[vid].e1)!=labs(points[vid].e5) &&
                       labs(points[vid].e1)!=labs(points[vid].e6) &&
                       labs(points[vid].e2)!=labs(points[vid].e3) &&
                       labs(points[vid].e2)!=labs(points[vid].e4) &&
                       labs(points[vid].e2)!=labs(points[vid].e5) &&
                       labs(points[vid].e2)!=labs(points[vid].e6) &&
                       labs(points[vid].e3)!=labs(points[vid].e4) &&
                       labs(points[vid].e3)!=labs(points[vid].e5) &&
                       labs(points[vid].e3)!=labs(points[vid].e6) &&
                       labs(points[vid].e4)!=labs(points[vid].e5) &&
                       labs(points[vid].e4)!=labs(points[vid].e6) &&
                       labs(points[vid].e5)!=labs(points[vid].e6));*/

                assert(labs(points[vid].e1)-1<getNumEdges());
                assert(labs(points[vid].e2)-1<getNumEdges());
                assert(labs(points[vid].e3)-1<getNumEdges());
                assert(labs(points[vid].e4)-1<getNumEdges());
                assert(labs(points[vid].e5)-1<getNumEdges());
                assert(labs(points[vid].e6)-1<getNumEdges());
            }
        }
    }
    checkInternalState();
}

bool DECMesh3D::checkInternalState()
{
    //Primal
    for(unsigned int i=0;i<getNumVoxels();i++)
    {
        Voxel3D v = voxels[i];
        assert(v.id<getNumVoxels());
        for(unsigned int j=0;j<6;j++)
        {
            assert(v.f[j]!=0);
        }

        for(unsigned int j=0;j<6;j+=2)
        {
            Face3D f1 = getFace(v.f[j+0]);
            Face3D f2 = getFace(v.f[j+1]);
            double s1 = getFaceSignum(v.f[j+0]);
            double s2 = getFaceSignum(v.f[j+1]);
            glm::dvec3 f1normal = glm::normalize(v.center-f1.center);
            glm::dvec3 f2normal = glm::normalize(v.center-f2.center);
            double dot1 = glm::dot(s1*f1normal,s2*f2normal);
            assert(dot1>0.0);
        }

        for(unsigned int j=0;j<2;j++)
        {
            Face3D f1 = getFace(v.f[j+0]);
            double s1 = getFaceSignum(v.f[j+0]);
            for(unsigned int k=0;k<2;k++)
            {
                Face3D f2 = getFace(v.f[j+2]);
                double s2 = getFaceSignum(v.f[j+2]);
                glm::dvec3 cross1 = glm::cross(s2*f2.normal,s1*f1.normal);
                assert(glm::dot(cross1,glm::dvec3(1.0,0.0,0.0))>0.0);
            }
        }

        for(unsigned int j=0;j<2;j++)
        {
            Face3D f1 = getFace(v.f[j+2]);
            double s1 = getFaceSignum(v.f[j+2]);
            for(unsigned int k=0;k<2;k++)
            {
                Face3D f2 = getFace(v.f[j+4]);
                double s2 = getFaceSignum(v.f[j+4]);
                glm::dvec3 cross1 = glm::cross(s2*f2.normal,s1*f1.normal);
                assert(glm::dot(cross1,glm::dvec3(0.0,0.0,1.0))>0.0);
            }
        }

        for(unsigned int j=0;j<2;j++)
        {
            Face3D f1 = getFace(v.f[j+4]);
            double s1 = getFaceSignum(v.f[j+4]);
            for(unsigned int k=0;k<2;k++)
            {
                Face3D f2 = getFace(v.f[j+0]);
                double s2 = getFaceSignum(v.f[j+0]);
                glm::dvec3 cross1 = glm::cross(s2*f2.normal,s1*f1.normal);
                assert(glm::dot(cross1,glm::dvec3(0.0,1.0,0.0))>0.0);
            }
        }
    }


    for(unsigned int i=0;i<getNumFaces();i++)
    {
        Face3D f = faces[i];
        assert(f.id<getNumFaces());
        for(unsigned int j=0;j<4;j++)
        {
            assert(f.e[j]!=0);
        }

        for(unsigned int j=0;j<4;j++)
        {
            Edge3D e1 = getEdge(f.e[(j+0)%4]);
            Edge3D e2 = getEdge(f.e[(j+1)%4]);

            Vertex3D e1v1 = getPoint(e1.v1);
            Vertex3D e1v2 = getPoint(e1.v2);
            Vertex3D e2v1 = getPoint(e2.v1);
            Vertex3D e2v2 = getPoint(e2.v2);

            glm::dvec3 evec1 = glm::dvec3(e1v2.center-e1v1.center);
            glm::dvec3 evec2 = glm::dvec3(e2v2.center-e2v1.center);

            double s1 = getFaceSignum(f.e[(j+0)%4]);
            double s2 = getFaceSignum(f.e[(j+1)%4]);

            glm::dvec3 cross1 = glm::normalize(glm::cross(s2*evec2,s1*evec1));
            assert(glm::dot(f.normal,cross1)<0.0);

        }
    }


    for(unsigned int i=0;i<getNumEdges();i++)
    {
        Edge3D e = edges[i];
        assert(e.id<getNumEdges());
        for(unsigned int j=0;j<2;j++)
        {
            assert(e.v[j]<getNumPoints());
        }
    }

    for (unsigned int i=0;i<getNumPoints();i++)
    {
        Vertex3D v = points[i];
        assert(v.id<getNumPoints());
    }

    //Dual
    for (unsigned int i=0;i<getNumPoints();i++)
    {
        Vertex3D v = points[i];
        assert(v.id<getNumPoints());

        unsigned int count=0;

        for(unsigned int j=0;j<2;j++)
        {
            if(v.e[j+0] == 0) continue;
            Edge3D f1 = getEdge(v.e[j+0]);
            double s1 = getFaceSignum(v.e[j+0]);
            Vertex3D e1v1 = getPoint(f1.v1);
            Vertex3D e1v2 = getPoint(f1.v2);
            glm::dvec3 f1normal = glm::normalize(e1v2.center-e1v1.center);

            for(unsigned int k=0;k<2;k++)
            {
                if(v.e[j+2] == 0) continue;
                Edge3D f2 = getEdge(v.e[j+2]);
                double s2 = getFaceSignum(v.e[j+2]);
                Vertex3D e2v1 = getPoint(f2.v1);
                Vertex3D e2v2 = getPoint(f2.v2);
                glm::dvec3 f2normal = glm::normalize(e2v2.center-e2v1.center);

                glm::dvec3 cross1 = glm::cross(s2*f2normal,s1*f1normal);
                assert(glm::dot(cross1,glm::dvec3(1.0,0.0,0.0))>0.0);
            }
        }

        for(unsigned int j=0;j<2;j++)
        {
            if(v.e[j+2] == 0) continue;
            Edge3D f1 = getEdge(v.e[j+2]);
            double s1 = getFaceSignum(v.e[j+2]);
            Vertex3D e1v1 = getPoint(f1.v1);
            Vertex3D e1v2 = getPoint(f1.v2);
            glm::dvec3 f1normal = glm::normalize(e1v2.center-e1v1.center);

            for(unsigned int k=0;k<2;k++)
            {
                if(v.e[j+4] == 0) continue;
                Edge3D f2 = getEdge(v.e[j+4]);
                double s2 = getFaceSignum(v.e[j+4]);
                Vertex3D e2v1 = getPoint(f2.v1);
                Vertex3D e2v2 = getPoint(f2.v2);
                glm::dvec3 f2normal = glm::normalize(e2v2.center-e2v1.center);

                glm::dvec3 cross1 = glm::cross(s2*f2normal,s1*f1normal);
                assert(glm::dot(cross1,glm::dvec3(0.0,0.0,1.0))>0.0);
            }
        }

        for(unsigned int j=0;j<2;j++)
        {
            if(v.e[j+4] == 0) continue;
            Edge3D f1 = getEdge(v.e[j+4]);
            double s1 = getFaceSignum(v.e[j+4]);
            Vertex3D e1v1 = getPoint(f1.v1);
            Vertex3D e1v2 = getPoint(f1.v2);
            glm::dvec3 f1normal = glm::normalize(e1v2.center-e1v1.center);

            for(unsigned int k=0;k<2;k++)
            {
                if(v.e[j+0] == 0) continue;
                Edge3D f2 = getEdge(v.e[j+0]);
                double s2 = getFaceSignum(v.e[j+0]);
                Vertex3D e2v1 = getPoint(f2.v1);
                Vertex3D e2v2 = getPoint(f2.v2);
                glm::dvec3 f2normal = glm::normalize(e2v2.center-e2v1.center);

                glm::dvec3 cross1 = glm::cross(s2*f2normal,s1*f1normal);
                assert(glm::dot(cross1,glm::dvec3(0.0,1.0,0.0))>0.0);
            }
        }

        for(unsigned int j=0;j<6;j+=2)
        {
            if(v.e[j+0]!=0 && v.e[j+1]!=0)
            {
                count++;
                Edge3D f1 = getEdge(v.e[j+0]);
                Edge3D f2 = getEdge(v.e[j+1]);
                double s1 = getFaceSignum(v.e[j+0]);
                double s2 = getFaceSignum(v.e[j+1]);
                glm::dvec3 f1normal = glm::normalize(v.center-f1.center);
                glm::dvec3 f2normal = glm::normalize(v.center-f2.center);
                double dot1 = glm::dot(s1*f1normal,s2*f2normal);
                assert(dot1>0.0);
            }
        }
    }

    //for(unsigned int i=0;i<getNumEdges();i++)
    for(unsigned int i=0;i<numZEdges;i++)
    {
        Edge3D f = edges[i];
        Vertex3D e1v1 = getPoint(f.v1);
        Vertex3D e1v2 = getPoint(f.v2);
        glm::dvec3 evec = glm::normalize(glm::dvec3(e1v2.center-e1v1.center));

        unsigned int count=0;
        if(f.f[0]!=0 && f.f[2]!=0)
        {
            count++;
            Face3D e1 = getFace(f.f[0]);
            Face3D e2 = getFace(f.f[2]);

            glm::dvec3 evec1 = glm::normalize(e1.normal);
            glm::dvec3 evec2 = glm::normalize(e2.normal);

            double s1 = getFaceSignum(f.f[0]);
            double s2 = getFaceSignum(f.f[2]);

            glm::dvec3 cross1 = glm::normalize(glm::cross(s2*evec2,s1*evec1));
            assert(abs(glm::dot(evec,cross1))>0.0);
        }

        if(f.f[1]!=0 && f.f[3]!=0)
        {
            count++;
            Face3D e1 = getFace(f.f[1]);
            Face3D e2 = getFace(f.f[3]);

            glm::dvec3 evec1 = glm::normalize(e1.normal);
            glm::dvec3 evec2 = glm::normalize(e2.normal);

            double s1 = getFaceSignum(f.f[1]);
            double s2 = getFaceSignum(f.f[3]);

            glm::dvec3 cross1 = glm::normalize(glm::cross(s2*evec2,s1*evec1));
            assert(abs(glm::dot(evec,cross1))>0.0);
        }

        //assert(count!=0);
    }

    return true;
}

void DECMesh3D::setPointInside(const Vertex3D& v)
{
    points[v.id].inside = v.inside;
}

void DECMesh3D::setEdgeInside(const Edge3D& e)
{
    unsigned int eid = labs(e.id)-1;
    assert(eid<getNumEdges());
    edges[eid].inside = e.inside;
    setPointInside(Vertex3D(edges[eid].v1,e.inside));
    setPointInside(Vertex3D(edges[eid].v2,e.inside));
}

void DECMesh3D::setFaceInside(const Face3D& f)
{
    unsigned int fid = labs(f.id)-1;
    assert(fid<getNumFaces());
    faces[fid].inside = f.inside;
    setEdgeInside(Edge3D(faces[fid].e1,f.inside));
    setEdgeInside(Edge3D(faces[fid].e2,f.inside));
    setEdgeInside(Edge3D(faces[fid].e3,f.inside));
    setEdgeInside(Edge3D(faces[fid].e4,f.inside));
}

void DECMesh3D::setVoxelInside(const Voxel3D& v)
{
    unsigned int vid = v.id;
    assert(vid<getNumVoxels());
    voxels[vid].inside = v.inside;
    setFaceInside(Face3D(voxels[vid].f1,v.inside));
    setFaceInside(Face3D(voxels[vid].f2,v.inside));
    setFaceInside(Face3D(voxels[vid].f3,v.inside));
    setFaceInside(Face3D(voxels[vid].f4,v.inside));
    setFaceInside(Face3D(voxels[vid].f5,v.inside));
    setFaceInside(Face3D(voxels[vid].f6,v.inside));
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

Voxel3D* DECMesh3D::getVoxels()
{
    return voxels.data();
}

Face3D* DECMesh3D::getFaces()
{
    return faces.data();
}

Edge3D* DECMesh3D::getEdges()
{
    return edges.data();
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

    return voxels[getVoxelIndex(xOfs,yOfs,zOfs)].inside==GridState::INSIDE;
}

Vertex3D DECMesh3D::getPoint(int id)
{
    return points[id];
}

Edge3D DECMesh3D::getEdge(int id)
{
    return edges[labs(id)-1];
}

Face3D DECMesh3D::getFace(int id)
{
    return faces[std::max(labs(id)-1,0l)];
}

Voxel3D DECMesh3D::getVoxel(int id)
{
    return voxels[id];
}

int DECMesh3D::getPointSignum(int id)
{
    return 1;
}

int DECMesh3D::getEdgeSignum(int id)
{
    return id>0?1:-1;
}

int DECMesh3D::getFaceSignum(int id)
{
    return id>0?1:-1;
}

int DECMesh3D::getVoxelSignum(int id)
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
