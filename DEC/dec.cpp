#include "dec.h"
#include <iostream>


Eigen::SparseMatrix<double> hodge2(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> h,b0,b1,b2;


    if(dual)
    {
        h.resize(mesh.getNumEdges(),mesh.getNumEdges());
        h.setIdentity();
        return h;

        b0 = derivative0(mesh);
        b1 = derivative1(mesh);
        b2 = derivative2(mesh);

        for(EdgeIterator eit=mesh.getEdgeIteratorBegin();eit<mesh.getEdgeIteratorEnd();eit++)
        {
            if(eit->inside==GridState::INSIDE)
            {
                unsigned int eid = mesh.getEdgeIndex(*eit);
                unsigned int nVoxels=0;
                for(Eigen::SparseMatrix<double>::InnerIterator it(b1,eid);it;++it)
                {
                    for(Eigen::SparseMatrix<double>::InnerIterator it2(b2,it.row());it2;++it2)
                    {
                        nVoxels++;
                    }
                }
                if(nVoxels==2)
                {
                    h.insert(eid,eid)=4.0/1.0;
                }
                else if(nVoxels==4)
                {
                    h.insert(eid,eid)=2.0/1.0;
                }
                else if(nVoxels==6)
                {
                    h.insert(eid,eid)=1.3333;
                }
                else if(nVoxels==8)
                {
                    h.insert(eid,eid)=1.0/1.0;
                }
            }
        }
    }
    else
    {
        h.resize(mesh.getNumFaces(),mesh.getNumFaces());
        h.setIdentity();
        return h;

        b0 = derivative0(mesh);
        b1 = derivative1(mesh);
        b2 = derivative2(mesh);

        for(FaceIterator fit=mesh.getFaceIteratorBegin();fit<mesh.getFaceIteratorEnd();fit++)
        {
            if(fit->inside==GridState::INSIDE)
            {
                unsigned int nVoxels=0;
                for(Eigen::SparseMatrix<double>::InnerIterator it(b1,mesh.getFaceIndex(*fit));it;++it)
                {
                    nVoxels++;
                }
                if(nVoxels==1)
                {
                    h.insert(mesh.getFaceIndex(*fit),mesh.getFaceIndex(*fit))=0.25/1.0;
                }
                else if(nVoxels==2)
                {
                    h.insert(mesh.getFaceIndex(*fit),mesh.getFaceIndex(*fit))=0.5/1.0;
                }
                else if(nVoxels==3)
                {
                    h.insert(mesh.getFaceIndex(*fit),mesh.getFaceIndex(*fit))=0.75/1.0;
                }
                else if(nVoxels==4)
                {
                    h.insert(mesh.getFaceIndex(*fit),mesh.getFaceIndex(*fit))=1.0/1.0;
                }
                else
                {
                    std::cout<<"Something mysterios"<<std::endl;
                }
            }
        }
    }

    return h;
}

Eigen::SparseMatrix<double> hodge1(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> h,d;
    d = derivative1(mesh);
    h.resize(mesh.getNumEdges(),mesh.getNumEdges());
    h.setIdentity();
    return h;
    for(EdgeIterator eit=mesh.getEdgeIteratorBegin();eit!=mesh.getEdgeIteratorEnd();eit++)
    {
        if(eit->inside==GridState::INSIDE)
        {
            unsigned int i=mesh.getEdgeIndex(*eit);
            unsigned int dualEdges=0;
            for(Eigen::SparseMatrix<double>::InnerIterator it(d,i);it;++it)
            {
                dualEdges++;
            }
            if(dual)
            {
                if(dualEdges==1)
                {
                    h.insert(i,i) = 2.0;
                }
                else if(dualEdges==2)
                {
                    h.insert(i,i) = 1.0;
                }
            }
            else
            {
                if(dualEdges==1)
                {
                    h.insert(i,i) = 0.5;
                }
                else if(dualEdges==2)
                {
                    h.insert(i,i) = 1.0;
                }
            }
        }
    }
    return h;
}

Eigen::SparseMatrix<double> hodge0(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> h;
    if(dual)
    {
        h.resize(mesh.getNumFaces(),mesh.getNumFaces());
    }
    else
    {
        h.resize(mesh.getNumPoints(),mesh.getNumPoints());
    }
    return h;
}

Eigen::SparseMatrix<double> derivative2(DECMesh3D& mesh,bool dual)
{
    if(dual)
    {
        return derivative0(mesh).transpose();
    }
    Eigen::SparseMatrix<double> d;
    d.resize(mesh.getNumFaces(),mesh.getNumVoxels());
    d.reserve(Eigen::VectorXi::Constant(mesh.getNumVoxels(),6));
    std::vector<Eigen::Triplet<double>> tripleList;
    tripleList.resize(mesh.getNumVoxels()*6);
    Voxel3D* voxels = mesh.getVoxels();
    #pragma omp parallel for
    for(int i=0;i<mesh.getNumVoxels();i++)
    {
        Voxel3D* voxel = &voxels[i];
        if(voxel->inside==GridState::INSIDE)
        {
            tripleList[i*6] = Eigen::Triplet<double>(labs(voxel->f1)-1,labs(voxel->id)-1, mesh.getFaceSignum(voxel->f1,voxel->v1,voxel->v2,voxel->v6,voxel->v5));
            tripleList[i*6+1] = Eigen::Triplet<double>(labs(voxel->f2)-1,labs(voxel->id)-1, mesh.getFaceSignum(voxel->f2,voxel->v8,voxel->v7,voxel->v3,voxel->v4));
            tripleList[i*6+2] = Eigen::Triplet<double>(labs(voxel->f3)-1,labs(voxel->id)-1, mesh.getFaceSignum(voxel->f3,voxel->v5,voxel->v6,voxel->v7,voxel->v8));
            tripleList[i*6+3] = Eigen::Triplet<double>(labs(voxel->f4)-1,labs(voxel->id)-1, mesh.getFaceSignum(voxel->f4,voxel->v4,voxel->v3,voxel->v2,voxel->v1));
            tripleList[i*6+4] = Eigen::Triplet<double>(labs(voxel->f5)-1,labs(voxel->id)-1, mesh.getFaceSignum(voxel->f5,voxel->v4,voxel->v1,voxel->v5,voxel->v8));
            tripleList[i*6+5] = Eigen::Triplet<double>(labs(voxel->f6)-1,labs(voxel->id)-1, mesh.getFaceSignum(voxel->f6,voxel->v2,voxel->v3,voxel->v7,voxel->v6));
        }
    }
    d.setFromTriplets(tripleList.begin(),tripleList.end());
    return d.transpose();
}

//TODO implement derivative
Eigen::SparseMatrix<double> derivative1(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> d;
    d.resize(mesh.getNumEdges(),mesh.getNumFaces());
    d.reserve(Eigen::VectorXi::Constant(mesh.getNumFaces(),4));
    std::vector<Eigen::Triplet<double>> tripleList;
    tripleList.resize(mesh.getNumFaces()*4);
    Face3D* faces = mesh.getFaces();
    #pragma omp parallel for
    for(int i=0;i<mesh.getNumFaces();i++)
    {
        Face3D* face = &faces[i];
        if(face->inside==GridState::INSIDE)
        {
            unsigned int v1 = face->v1;
            unsigned int v2 = face->v2;
            unsigned int v3 = face->v3;
            unsigned int v4 = face->v4;

            tripleList[i*4] = Eigen::Triplet<double>(face->e1,mesh.getFaceIndex(*face),mesh.getEdgeSignum(face->e1,v1,v2));
            tripleList[i*4+1] = Eigen::Triplet<double>(face->e2,mesh.getFaceIndex(*face),mesh.getEdgeSignum(face->e2,v2,v3));
            tripleList[i*4+2] = Eigen::Triplet<double>(face->e3,mesh.getFaceIndex(*face),mesh.getEdgeSignum(face->e3,v3,v4));
            tripleList[i*4+3] = Eigen::Triplet<double>(face->e4,mesh.getFaceIndex(*face),mesh.getEdgeSignum(face->e4,v4,v1));
        }
    }
    d.setFromTriplets(tripleList.begin(),tripleList.end());
    if(dual)
    {
        return d;
    }
    else
    {
        return d.transpose();
    }

}

Eigen::SparseMatrix<double> derivative0(DECMesh3D& mesh,bool dual)
{
    if(dual)
    {
        return derivative2(mesh).transpose();
    }
    Eigen::SparseMatrix<double> d;
    d.resize(mesh.getNumPoints(),mesh.getNumEdges());
    d.reserve(Eigen::VectorXi::Constant(mesh.getNumEdges(),2));
    for(EdgeIterator it = mesh.getEdgeIteratorBegin();it!=mesh.getEdgeIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            unsigned int eid = mesh.getEdgeIndex(*it);
            unsigned int v1 = it->v1;
            unsigned int v2 = it->v2;
            d.insert(v1,eid) = -1.0;
            d.insert(v2,eid) = 1.0;
        }
    }
    return d.transpose();
}
