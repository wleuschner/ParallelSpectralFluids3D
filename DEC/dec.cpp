#include "dec.h"
#include <iostream>


Eigen::SparseMatrix<double> hodge2(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> h;
    std::vector<Eigen::Triplet<double>> tripletList;

    if(dual)
    {
        h.resize(mesh.getNumEdges(),mesh.getNumEdges());

        tripletList.resize(mesh.getNumEdges());

        for(PointIterator vit=mesh.getPointIteratorBegin();vit<mesh.getPointIteratorEnd();vit++)
        {
            if(vit->inside==GridState::INSIDE)
            {
                for(unsigned int i=0;i<6;i++)
                {
                    if(vit->e[i]==0) continue;
                    Edge3D f = mesh.getEdge(vit->e[i]);
                    if(f.inside!=GridState::INSIDE) continue;
                    glm::dvec3 v1 = mesh.getPoint(f.v1).center;
                    glm::dvec3 v2 = mesh.getPoint(f.v2).center;
                    glm::dvec3 normal1= glm::normalize(v2-v1);
                    double areaPrim = 0.0;
                    double areaPrim3 = 0.0;
                    for(unsigned int j=0;j<4;j++)
                    {
                        Face3D e1 = mesh.getFace(f.f[j]);
                        if(e1.inside!=GridState::INSIDE) continue;
                        for(unsigned int k=0;k<4;k++)
                        {
                            Face3D e2 = mesh.getFace(f.f[k]);
                            if(e2.inside!=GridState::INSIDE) continue;
                            if(e2.id==e1.id) continue;

                            double side1 = glm::length(f.center-e1.center);
                            double side2 = glm::length(f.center-e2.center);

                            glm::dvec3 cross1 = glm::cross(e1.normal,e2.normal);
                            areaPrim += side1*side2*abs(glm::dot(normal1,cross1))*(f.f[j]!=0)*(f.f[k]!=0);
                        }
                    }
                    areaPrim/=2;

                    /*if(f.f1!=0 && f.f2!=0)
                    {
                        Face3D e1 = mesh.getFace(f.f1);
                        Face3D e2 = mesh.getFace(f.f2);

                        {
                            double side1 = glm::length(f.center-e1.center);
                            double side2 = glm::length(f.center-e2.center);
                            double area = side1*side2;
                            areaPrim += area;
                            assert(abs(double(mesh.resolution)/2-side1)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)/2-side2)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                        }
                    }
                    if(f.f2!=0 && f.f3!=0)
                    {
                        Face3D e2 = mesh.getFace(f.f2);
                        Face3D e3 = mesh.getFace(f.f3);
                        {
                            double side1 = glm::length(f.center-e2.center);
                            double side2 = glm::length(f.center-e3.center);
                            double area = side1*side2;
                            areaPrim += area;
                            assert(abs(double(mesh.resolution)/2-side1)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)/2-side2)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                        }
                    }
                    if(f.f3!=0 && f.f4!=0)
                    {
                        Face3D e3 = mesh.getFace(f.f3);
                        Face3D e4 = mesh.getFace(f.f4);

                        {
                            double side1 = glm::length(f.center-e3.center);
                            double side2 = glm::length(f.center-e4.center);
                            double area = side1*side2;
                            areaPrim += area;
                            assert(abs(double(mesh.resolution)/2-side1)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)/2-side2)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                        }
                    }
                    if(f.f4!=0 && f.f1!=0)
                    {
                        Face3D e1 = mesh.getFace(f.f1);
                        Face3D e4 = mesh.getFace(f.f4);

                        {
                            double side1 = glm::length(f.center-e4.center);
                            double side2 = glm::length(f.center-e1.center);
                            double area = side1*side2;
                            areaPrim += area;
                            assert(abs(double(mesh.resolution)/2-side1)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)/2-side2)<=std::numeric_limits<float>::epsilon());
                            assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                        }
                    }*/

                    double areaPrim2 = 0;
                    if(f.dualCount==2)
                    {
                        areaPrim2 = 1*(double(mesh.resolution)/2*double(mesh.resolution)/2);
                    }
                    else if(f.dualCount==3)
                    {
                        areaPrim2 = 2*(double(mesh.resolution)/2*double(mesh.resolution)/2);
                    }
                    else if(f.dualCount==4)
                    {
                        areaPrim2 = 4*(double(mesh.resolution)/2*double(mesh.resolution)/2);
                    }
                    assert(abs(areaPrim-areaPrim2)<std::numeric_limits<float>::epsilon());
                    assert(areaPrim!=0.0);
                    if(areaPrim!=0.0)
                    {
                        const_cast<int&>(tripletList[labs(vit->e[i])-1].row())=labs(vit->e[i])-1;
                        const_cast<int&>(tripletList[labs(vit->e[i])-1].col())=labs(vit->e[i])-1;
                        const_cast<double&>(tripletList[labs(vit->e[i])-1].value())=1.0/(glm::length((vit->center-f.center))/areaPrim);
                    }
                    else {
                       std::cout<<f.id<<" "<<f.f1<<" "<<f.f2<<" "<<f.f3<<" "<<f.f4<<std::endl;
                       //exit(-1);
                    }
                }
            }
        }
        h.setFromTriplets(tripletList.begin(),tripletList.end());
    }
    else
    {
        h.resize(mesh.getNumFaces(),mesh.getNumFaces());
        tripletList.resize(mesh.getNumFaces());

        for(VoxelIterator vit=mesh.getVoxelIteratorBegin();vit<mesh.getVoxelIteratorEnd();vit++)
        {
                for(unsigned int i=0;i<6;i++)
                {
                    Face3D f = mesh.getFace(vit->f[i]);
                    double areaPrim = 0.0;
                    Edge3D e1 = mesh.getEdge(f.e1);
                    Edge3D e2 = mesh.getEdge(f.e2);
                    Edge3D e3 = mesh.getEdge(f.e3);
                    Edge3D e4 = mesh.getEdge(f.e4);

                    assert(e1.id!=e2.id &&
                           e1.id!=e3.id &&
                           e1.id!=e4.id &&
                           e2.id!=e3.id &&
                           e2.id!=e4.id &&
                           e3.id!=e4.id);

                    {
                        double area = glm::length(f.center-e1.center)*glm::length(f.center-e2.center);
                        areaPrim += area;
                        assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                    }
                    {
                        double area = glm::length(f.center-e2.center)*glm::length(f.center-e3.center);
                        areaPrim += area;
                        assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                    }
                    {
                        double area = glm::length(f.center-e3.center)*glm::length(f.center-e4.center);
                        areaPrim += area;
                        assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                    }
                    {
                        double area = glm::length(f.center-e4.center)*glm::length(f.center-e1.center);
                        areaPrim += area;
                        assert(abs(double(mesh.resolution)*double(mesh.resolution)/4-area)<=std::numeric_limits<float>::epsilon());
                    }

                    assert(areaPrim!=0.0);
                    const_cast<int&>(tripletList[labs(vit->f[i])-1].row())=labs(vit->f[i])-1;
                    const_cast<int&>(tripletList[labs(vit->f[i])-1].col())=labs(vit->f[i])-1;
                    const_cast<double&>(tripletList[labs(vit->f[i])-1].value())=((glm::length(vit->center-f.center))/areaPrim);
                }
        }
        h.setFromTriplets(tripletList.begin(),tripletList.end());
    }

    h.setIdentity();
    return h;
}

Eigen::SparseMatrix<double> hodge1(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> h,d;
    d = derivative1(mesh);
    h.resize(mesh.getNumEdges(),mesh.getNumEdges());

    for(FaceIterator fit=mesh.getFaceIteratorBegin();fit!=mesh.getFaceIteratorEnd();fit++)
    {
        if(fit->inside==GridState::INSIDE)
        {
            Edge3D e1 = mesh.getEdge(fit->e1);
            Edge3D e2 = mesh.getEdge(fit->e2);
            Edge3D e3 = mesh.getEdge(fit->e3);
            Edge3D e4 = mesh.getEdge(fit->e4);

            std::vector<Eigen::Triplet<double>> tripletList;
            tripletList.resize(mesh.getNumEdges());

            if(dual)
            {
                if(fit->e4!=0)
                {
                    const_cast<int&>(tripletList[labs(fit->e1)-1].row())=labs(fit->e1)-1;
                    const_cast<int&>(tripletList[labs(fit->e1)-1].col())=labs(fit->e1)-1;
                    const_cast<double&>(tripletList[labs(fit->e1)-1].value())+=glm::length(fit->center-e1.center);
                }

                if(fit->e4!=0)
                {
                    const_cast<int&>(tripletList[labs(fit->e2)-1].row())=labs(fit->e2)-1;
                    const_cast<int&>(tripletList[labs(fit->e2)-1].col())=labs(fit->e2)-1;
                    const_cast<double&>(tripletList[labs(fit->e2)-1].value())+=glm::length(fit->center-e2.center);
                }

                if(fit->e4!=0)
                {
                    const_cast<int&>(tripletList[labs(fit->e3)-1].row())=labs(fit->e3)-1;
                    const_cast<int&>(tripletList[labs(fit->e3)-1].col())=labs(fit->e3)-1;
                    const_cast<double&>(tripletList[labs(fit->e3)-1].value())+=glm::length(fit->center-e3.center);
                }

                if(fit->e4!=0)
                {
                    const_cast<int&>(tripletList[labs(fit->e4)-1].row())=labs(fit->e4)-1;
                    const_cast<int&>(tripletList[labs(fit->e4)-1].col())=labs(fit->e4)-1;
                    const_cast<double&>(tripletList[labs(fit->e4)-1].value())+=glm::length(fit->center-e4.center);
                }
            }
            else
            {
                const_cast<int&>(tripletList[labs(fit->e1)-1].row())=labs(fit->e1)-1;
                const_cast<int&>(tripletList[labs(fit->e1)-1].col())=labs(fit->e1)-1;
                const_cast<double&>(tripletList[labs(fit->e1)-1].value())+=1.0/glm::length(fit->center-e1.center);

                const_cast<int&>(tripletList[labs(fit->e2)-1].row())=labs(fit->e2)-1;
                const_cast<int&>(tripletList[labs(fit->e2)-1].col())=labs(fit->e2)-1;
                const_cast<double&>(tripletList[labs(fit->e2)-1].value())+=1.0/glm::length(fit->center-e2.center);

                const_cast<int&>(tripletList[labs(fit->e3)-1].row())=labs(fit->e3)-1;
                const_cast<int&>(tripletList[labs(fit->e3)-1].col())=labs(fit->e3)-1;
                const_cast<double&>(tripletList[labs(fit->e3)-1].value())+=1.0/glm::length(fit->center-e3.center);

                const_cast<int&>(tripletList[labs(fit->e4)-1].row())=labs(fit->e4)-1;
                const_cast<int&>(tripletList[labs(fit->e4)-1].col())=labs(fit->e4)-1;
                const_cast<double&>(tripletList[labs(fit->e4)-1].value())+=1.0/glm::length(fit->center-e4.center);
            }

            /*
            unsigned int i=mesh.getEdgeIndex(*it);
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
            */
        }
    }

    /*
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
    }*/
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
            tripleList[i*6] = Eigen::Triplet<double>(labs(voxel->f1)-1,voxel->id, voxel->f1>0?1:-1);
            tripleList[i*6+1] = Eigen::Triplet<double>(labs(voxel->f2)-1,voxel->id, voxel->f2>0?1:-1);
            tripleList[i*6+2] = Eigen::Triplet<double>(labs(voxel->f3)-1,voxel->id, voxel->f3>0?1:-1);
            tripleList[i*6+3] = Eigen::Triplet<double>(labs(voxel->f4)-1,voxel->id, voxel->f4>0?1:-1);
            tripleList[i*6+4] = Eigen::Triplet<double>(labs(voxel->f5)-1,voxel->id, voxel->f5>0?1:-1);
            tripleList[i*6+5] = Eigen::Triplet<double>(labs(voxel->f6)-1,voxel->id, voxel->f6>0?1:-1);
        }
    }
    d.setFromTriplets(tripleList.begin(),tripleList.end());
    return d.transpose();
}

//TODO implement derivative
Eigen::SparseMatrix<double> derivative1(DECMesh3D& mesh,bool dual)
{
    Eigen::SparseMatrix<double> d;

    //if(dual)
    {
        /*d.resize(mesh.getNumEdges(),mesh.getNumFaces());
        d.reserve(Eigen::VectorXi::Constant(mesh.getNumFaces(),4));
        std::vector<Eigen::Triplet<double>> tripleList;
        tripleList.resize(mesh.getNumFaces()*4);
        Face3D* faces = mesh.getFaces();
        //#pragma omp parallel for
        for(int i=0;i<mesh.getNumFaces();i++)
        {
            Face3D* face = &faces[i];
            if(face->inside == GridState::INSIDE)
            {
                tripleList[i*4] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e1),mesh.getFaceIndex(*face),(face->e1>=0?1:-1));
                tripleList[i*4+1] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e2),mesh.getFaceIndex(*face),(face->e2>=0?1:-1));
                tripleList[i*4+2] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e3),mesh.getFaceIndex(*face),(face->e3>=0?1:-1));
                tripleList[i*4+3] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e4),mesh.getFaceIndex(*face),(face->e4>=0?1:-1));
            }
        }
        d.setFromTriplets(tripleList.begin(),tripleList.end());*/
    }
    //else
    {
        d.resize(mesh.getNumEdges(),mesh.getNumFaces());
        d.reserve(Eigen::VectorXi::Constant(mesh.getNumFaces(),4));
        std::vector<Eigen::Triplet<double>> tripleList;
        tripleList.resize(mesh.getNumFaces()*4);
        Face3D* faces = mesh.getFaces();
        //#pragma omp parallel for
        for(int i=0;i<mesh.getNumFaces();i++)
        {
            Face3D* face = &faces[i];
            if(face->inside == GridState::INSIDE)
            {
                tripleList[i*4] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e1),mesh.getFaceIndex(*face),(face->e1>=0?1:-1));
                tripleList[i*4+1] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e2),mesh.getFaceIndex(*face),(face->e2>=0?1:-1));
                tripleList[i*4+2] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e3),mesh.getFaceIndex(*face),(face->e3>=0?1:-1));
                tripleList[i*4+3] = Eigen::Triplet<double>(mesh.signedIdToIndex(face->e4),mesh.getFaceIndex(*face),(face->e4>=0?1:-1));
            }
        }

/*
        Voxel3D* voxels = mesh.getVoxels();
        for(unsigned int k=0;k<mesh.getNumVoxels();k++)
        {
            Voxel3D* voxel = &voxels[k];
            for(int i=0;i<6;i++)
            {
                Face3D face2 = mesh.getFace(voxel->f[i]);
                Face3D* face = &face2;
                double fs = mesh.getFaceSignum(voxel->f[i]);
                if(face->inside == GridState::INSIDE)
                {
                    const_cast<double&>(tripleList[face->id*4].value()) *= fs;
                    const_cast<double&>(tripleList[face->id*4+1].value()) *= fs;
                    const_cast<double&>(tripleList[face->id*4+2].value()) *= fs;
                    const_cast<double&>(tripleList[face->id*4+3].value()) *= fs;
                }
            }
        }*/
        d.setFromTriplets(tripleList.begin(),tripleList.end());
    }
    if(dual)
    {
        return d;
    }
    else
    {
    }
    return d.transpose();
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
