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
                unsigned int nVoxels=0;
                for(Eigen::SparseMatrix<double>::InnerIterator it(b1,eit->id);it;++it)
                {
                    for(Eigen::SparseMatrix<double>::InnerIterator it2(b2,it.row());it2;++it2)
                    {
                        nVoxels++;
                    }
                }
                if(nVoxels==2)
                {
                    h.insert(eit->id,eit->id)=1.0/1.0;
                }
                else if(nVoxels==4)
                {
                    h.insert(eit->id,eit->id)=1.0/1.0;
                }
                else if(nVoxels==6)
                {
                    h.insert(eit->id,eit->id)=1.0/1.0;
                }
                else if(nVoxels==8)
                {
                    h.insert(eit->id,eit->id)=1.0/1.0;
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
                for(Eigen::SparseMatrix<double>::InnerIterator it(b2,fit->id);it;++it)
                {
                    nVoxels++;
                }
                if(nVoxels==1)
                {
                    h.insert(fit->id,fit->id)=1.0/1.0;
                }
                else if(nVoxels==2)
                {
                    h.insert(fit->id,fit->id)=1.0/1.0;
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
    for(EdgeIterator eit=mesh.getEdgeIteratorBegin();eit!=mesh.getEdgeIteratorEnd();eit++)
    {
        if(eit->inside==GridState::INSIDE)
        {
            unsigned int i=eit->id;
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
    for(VoxelIterator it = mesh.getVoxelIteratorBegin();it!=mesh.getVoxelIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            d.insert(it->f1,it->id) = mesh.getFaceSignum(it->f1,it->v1,it->v2,it->v6,it->v5);
            d.insert(it->f2,it->id) = mesh.getFaceSignum(it->f2,it->v8,it->v7,it->v3,it->v4);
            d.insert(it->f3,it->id) = mesh.getFaceSignum(it->f3,it->v5,it->v6,it->v7,it->v8);
            d.insert(it->f4,it->id) = mesh.getFaceSignum(it->f4,it->v4,it->v3,it->v2,it->v1);
            d.insert(it->f5,it->id) = mesh.getFaceSignum(it->f5,it->v4,it->v1,it->v5,it->v8);
            d.insert(it->f6,it->id) = mesh.getFaceSignum(it->f6,it->v2,it->v3,it->v7,it->v6);
        }
    }
    return d.transpose();
}

//TODO implement derivative
Eigen::SparseMatrix<double> derivative1(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> d;
    d.resize(mesh.getNumEdges(),mesh.getNumFaces());
    for(FaceIterator it = mesh.getFaceIteratorBegin();it!=mesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            unsigned int v1 = it->v1;
            unsigned int v2 = it->v2;
            unsigned int v3 = it->v3;
            unsigned int v4 = it->v4;

            d.insert(it->e1,it->id) = mesh.getEdgeSignum(it->e1,v1,v2);
            d.insert(it->e2,it->id) = mesh.getEdgeSignum(it->e2,v2,v3);
            d.insert(it->e3,it->id) = mesh.getEdgeSignum(it->e3,v3,v4);
            d.insert(it->e4,it->id) = mesh.getEdgeSignum(it->e4,v4,v1);

        }
    }
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
            unsigned int v1 = it->v1;
            unsigned int v2 = it->v2;
            d.insert(v1,it->id) = -1.0;
            d.insert(v2,it->id) = 1.0;
        }
    }
    return d.transpose();
}
