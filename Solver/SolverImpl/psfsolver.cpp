#include"psfsolver.h"
#include <glm/gtc/matrix_transform.hpp>
#include<GL/glew.h>
#include<iostream>
#include<set>
#include<map>
#include"../../Spectra/MatOp/SparseSymShiftSolve.h"
#include"../../Spectra/SymEigsShiftSolver.h"
#include"../../Graphics/TextureArray/TextureArray.h"
#include"../../Graphics/FrameBufferObject/FrameBufferObject.h"
#include"../../DEC/dec.h"

PSFSolver::PSFSolver() : AbstractSolver()
{
}

void PSFSolver::integrate()
{
    double e1 = 0.0;
    double e2 = 0.0;

    e1 = basisCoeff.dot(basisCoeff);

    Eigen::VectorXd vel(nEigenFunctions);
    #pragma omp parallel for
    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        vel(k) = (basisCoeff.transpose()*advection[k]*basisCoeff);
    }
    basisCoeff += timeStep*vel;

    e2 = basisCoeff.dot(basisCoeff);

    basisCoeff *= std::sqrt(e1/e2);

    #pragma omp parallel for
    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        basisCoeff(k) *= std::exp(-viscosity*eigenValues(k)*(timeStep));
    }
    if(gravityActive)
    {
        basisCoeff += timeStep*gravity;
    }


    velocityField = velBasisField*basisCoeff;

    glm::uvec3 dims = decMesh.getDimensions();
    std::vector<Vertex>& vertex = mesh->getVertices();
    #pragma omp parallel for
    for(std::vector<Particle>::iterator it=particles.begin();it<particles.end();it++)
    {
        if(simTime<it->lifeTime)
        {
            unsigned int yOfs = static_cast<unsigned int>((-mesh->getAABB().min.y+it->position.y)/resolution)+1;
            unsigned int xOfs = static_cast<unsigned int>((-mesh->getAABB().min.x+it->position.x)/resolution)+1;
            unsigned int zOfs = static_cast<unsigned int>((-mesh->getAABB().min.z+it->position.z)/resolution)+1;

            unsigned int zOfsMinus = static_cast<unsigned int>((-mesh->getAABB().min.z+it->position.z-resolution/2)/resolution)+1;
            unsigned int zOfsPlus = static_cast<unsigned int>((-mesh->getAABB().min.z+it->position.z+resolution/2)/resolution)+1;
            unsigned int yOfsMinus = static_cast<unsigned int>((-mesh->getAABB().min.y+it->position.y-resolution/2)/resolution)+1;
            unsigned int yOfsPlus = static_cast<unsigned int>((-mesh->getAABB().min.y+it->position.y+resolution/2)/resolution)+1;
            unsigned int xOfsMinus = static_cast<unsigned int>((-mesh->getAABB().min.x+it->position.x-resolution/2)/resolution)+1;
            unsigned int xOfsPlus = static_cast<unsigned int>((-mesh->getAABB().min.x+it->position.x+resolution/2)/resolution)+1;

            Voxel3D v1x = decMesh.getVoxel(zOfsMinus*(dims.y*dims.x)+yOfsMinus*dims.x+xOfs);
            Voxel3D v2x = decMesh.getVoxel(zOfsMinus*(dims.y*dims.x)+yOfsPlus*dims.x+xOfs);
            Voxel3D v3x = decMesh.getVoxel(zOfsPlus*(dims.y*dims.x)+yOfsMinus*dims.x+xOfs);
            Voxel3D v4x = decMesh.getVoxel(zOfsPlus*(dims.y*dims.x)+yOfsPlus*dims.x+xOfs);

            Voxel3D v1y = decMesh.getVoxel(zOfsMinus*(dims.y*dims.x)+yOfs*dims.x+xOfsMinus);
            Voxel3D v2y = decMesh.getVoxel(zOfsMinus*(dims.y*dims.x)+yOfs*dims.x+xOfsPlus);
            Voxel3D v3y = decMesh.getVoxel(zOfsPlus*(dims.y*dims.x)+yOfs*dims.x+xOfsMinus);
            Voxel3D v4y = decMesh.getVoxel(zOfsPlus*(dims.y*dims.x)+yOfs*dims.x+xOfsPlus);

            Voxel3D v1z = decMesh.getVoxel(zOfs*(dims.y*dims.x)+yOfsMinus*dims.x+xOfsMinus);
            Voxel3D v2z = decMesh.getVoxel(zOfs*(dims.y*dims.x)+yOfsMinus*dims.x+xOfsPlus);
            Voxel3D v3z = decMesh.getVoxel(zOfs*(dims.y*dims.x)+yOfsPlus*dims.x+xOfsMinus);
            Voxel3D v4z = decMesh.getVoxel(zOfs*(dims.y*dims.x)+yOfsPlus*dims.x+xOfsPlus);

            glm::vec3 cf1x = 0.25f*(vertex[v1x.v5].pos+vertex[v1x.v1].pos+vertex[v1x.v5].pos+vertex[v1x.v8].pos);
            glm::vec3 cf1y = 0.25f*(vertex[v1y.v7].pos+vertex[v1y.v8].pos+vertex[v1y.v5].pos+vertex[v1y.v6].pos);
            glm::vec3 cf1z = 0.25f*(vertex[v1z.v1].pos+vertex[v1z.v2].pos+vertex[v1z.v6].pos+vertex[v1z.v5].pos);

            double vel1x,vel2x,vel3x,vel4x,vel5x,vel6x,vel7x,vel8x;
            vel1x = decMesh.getFaceSignum(v1x.f5,v1x.v4,v1x.v1,v1x.v5,v1x.v8)*velocityField(v1x.f5);
            vel3x = -decMesh.getFaceSignum(v1x.f6,v1x.v2,v1x.v3,v1x.v7,v1x.v6)*velocityField(v1x.f6);
            vel2x = decMesh.getFaceSignum(v2x.f5,v2x.v4,v2x.v1,v2x.v5,v2x.v8)*velocityField(v2x.f5);
            vel4x = -decMesh.getFaceSignum(v2x.f6,v2x.v2,v2x.v3,v2x.v7,v2x.v6)*velocityField(v2x.f6);
            vel5x = decMesh.getFaceSignum(v3x.f5,v3x.v4,v3x.v1,v3x.v5,v3x.v8)*velocityField(v3x.f5);
            vel7x = -decMesh.getFaceSignum(v3x.f6,v3x.v2,v3x.v3,v3x.v7,v3x.v6)*velocityField(v3x.f6);
            vel6x = decMesh.getFaceSignum(v4x.f5,v4x.v4,v4x.v1,v4x.v5,v4x.v8)*velocityField(v4x.f5);
            vel8x = -decMesh.getFaceSignum(v4x.f6,v4x.v2,v4x.v3,v4x.v7,v4x.v6)*velocityField(v4x.f6);

            double vel1y,vel2y,vel3y,vel4y,vel5y,vel6y,vel7y,vel8y;
            vel1y = decMesh.getFaceSignum(v1y.f3,v1y.v7,v1y.v8,v1y.v5,v1y.v6)*velocityField(v1y.f3);
            vel2y = -decMesh.getFaceSignum(v1y.f4,v1y.v4,v1y.v3,v1y.v2,v1y.v1)*velocityField(v1y.f4);
            vel3y = decMesh.getFaceSignum(v2y.f3,v2y.v7,v2y.v8,v2y.v5,v2y.v6)*velocityField(v2y.f3);
            vel4y = -decMesh.getFaceSignum(v2y.f4,v2y.v4,v2y.v3,v2y.v2,v2y.v1)*velocityField(v2y.f4);
            vel5y = decMesh.getFaceSignum(v3y.f3,v3y.v7,v3y.v8,v3y.v5,v3y.v6)*velocityField(v3y.f3);
            vel6y = -decMesh.getFaceSignum(v3y.f4,v3y.v4,v3y.v3,v3y.v2,v3y.v1)*velocityField(v3y.f4);
            vel7y = decMesh.getFaceSignum(v4y.f3,v4y.v7,v4y.v8,v4y.v5,v4y.v6)*velocityField(v4y.f3);
            vel8y = -decMesh.getFaceSignum(v4y.f4,v4y.v4,v4y.v3,v4y.v2,v4y.v1)*velocityField(v4y.f4);

            double vel1z,vel2z,vel3z,vel4z,vel5z,vel6z,vel7z,vel8z;
            vel1z = decMesh.getFaceSignum(v1z.f1,v1z.v1,v1z.v2,v1z.v6,v1z.v5)*velocityField(v1z.f1);
            vel2z = decMesh.getFaceSignum(v3z.f1,v3z.v1,v3z.v2,v3z.v6,v3z.v5)*velocityField(v3z.f1);
            vel3z = decMesh.getFaceSignum(v2z.f1,v2z.v1,v2z.v2,v2z.v6,v2z.v5)*velocityField(v2z.f1);
            vel4z = decMesh.getFaceSignum(v4z.f1,v4z.v1,v4z.v2,v4z.v6,v4z.v5)*velocityField(v4z.f1);
            vel5z = -decMesh.getFaceSignum(v1z.f2,v1z.v3,v1z.v4,v1z.v8,v1z.v7)*velocityField(v1z.f2);
            vel6z = -decMesh.getFaceSignum(v3z.f2,v3z.v3,v3z.v4,v3z.v8,v3z.v7)*velocityField(v3z.f2);
            vel7z = -decMesh.getFaceSignum(v2z.f2,v2z.v3,v2z.v4,v2z.v8,v2z.v7)*velocityField(v2z.f2);
            vel8z = -decMesh.getFaceSignum(v4z.f2,v4z.v3,v4z.v4,v4z.v8,v4z.v7)*velocityField(v4z.f2);

            glm::dvec3 particleNormalizedX;
            glm::dvec3 particleNormalizedY;
            glm::dvec3 particleNormalizedZ;
            glm::dvec3 vel;

            particleNormalizedX = (1.0/resolution)*((it->position)-glm::dvec3(cf1x));
            glm::dvec4 xVelXInterp = glm::mix(glm::dvec4(vel1x,vel3x,vel5x,vel7x),glm::dvec4(vel2x,vel4x,vel6x,vel8x),particleNormalizedX.y);
            glm::dvec2 xVelYInterp = glm::mix(glm::dvec2(xVelXInterp.x,xVelXInterp.z),glm::dvec2(xVelXInterp.y,xVelXInterp.w),particleNormalizedX.x);
            vel.x = glm::mix(xVelYInterp.x,xVelYInterp.y,particleNormalizedX.z);

            particleNormalizedY = (1.0/resolution)*((it->position)-glm::dvec3(cf1y));
            glm::dvec4 yVelXInterp = glm::mix(glm::dvec4(vel1y,vel3y,vel5y,vel7y),glm::dvec4(vel2y,vel4y,vel6y,vel8y),particleNormalizedY.y);
            glm::dvec2 yVelYInterp = glm::mix(glm::dvec2(yVelXInterp.x,yVelXInterp.z),glm::dvec2(yVelXInterp.y,yVelXInterp.w),particleNormalizedY.x);
            vel.y = glm::mix(yVelYInterp.x,yVelYInterp.y,particleNormalizedY.z);

            particleNormalizedZ = (1.0/resolution)*((it->position)-glm::dvec3(cf1z));
            glm::dvec4 zVelXInterp = glm::mix(glm::dvec4(vel1z,vel3z,vel5z,vel7z),glm::dvec4(vel2z,vel4z,vel6z,vel8z),particleNormalizedZ.y);
            glm::dvec2 zVelYInterp = glm::mix(glm::dvec2(zVelXInterp.x,zVelXInterp.z),glm::dvec2(zVelXInterp.y,zVelXInterp.w),particleNormalizedZ.x);
            vel.z = glm::mix(zVelYInterp.x,zVelYInterp.y,particleNormalizedZ.z);

            //it->position = (it->position)+timeStep*vel;
            it->position = glm::clamp(it->position+timeStep*vel,glm::dvec3(mesh->getAABB().min),glm::dvec3(mesh->getAABB().max));
        }
        else
        {
            it->lifeTime = 0.0;
        }
    }
    simTime+=1.0/60.0;
}

void PSFSolver::buildLaplace()
{
    Eigen::SparseMatrix<double> mat = -1.0*derivative1(decMesh,false)*hodge2(decMesh,true)*derivative1(decMesh,true)*hodge2(decMesh,false);
    Eigen::SparseMatrix<double> bound = derivative2(decMesh);
    curl = derivative1(decMesh,true)*hodge2(decMesh,false);
    for(int k=0;k<bound.outerSize();k++)
    {
        unsigned int nVoxel=0;
        for(Eigen::SparseMatrix<double>::InnerIterator it(bound,k);it;++it)
        {
            nVoxel++;
        }
        if(nVoxel!=2)
        {
            //mat.prune([k](int i,int j,double v){return !(i==k||j==k);});
        }
    }
    mat.pruned();

    eigenValues.resize(nEigenFunctions);
    velBasisField.resize(decMesh.getNumFaces(),nEigenFunctions);
    bool decompositionDone=false;
    unsigned int foundEigenValues=0;
    double omega = 1e-6;
    while(foundEigenValues<nEigenFunctions)
    {
        try
        {
            Spectra::SparseSymShiftSolve<double> op(mat);
            Spectra::SymEigsShiftSolver<double,Spectra::WHICH_LM,Spectra::SparseSymShiftSolve<double>> solver(&op,nEigenFunctions,2*nEigenFunctions,omega);
            solver.init();

            int nconv = solver.compute(1000,1e-6,Spectra::WHICH_SM);
            if(solver.info()==Spectra::SUCCESSFUL)
            {

                Eigen::VectorXd tempEigenValues = solver.eigenvalues().real();
                Eigen::MatrixXd tempEigenVectors = solver.eigenvectors().real();
                for(unsigned int i=0;i<nconv && foundEigenValues<nEigenFunctions;i++)
                {
                    if(std::abs(tempEigenValues(i))>1e-10)
                    {
                        for(unsigned int j=0;j<foundEigenValues;j++)
                        {
                            if((velBasisField.col(j).dot(tempEigenVectors.col(i)))>1.0-std::numeric_limits<double>::epsilon())
                            {
                                std::cout<<"Already Inside"<<std::endl;
                                continue;
                            }
                        }
                        std::cout<<"ONE GARBAGE "<<tempEigenValues(i)<<std::endl;
                        omega = tempEigenValues(i);
                        eigenValues(foundEigenValues) = tempEigenValues(i);
                        velBasisField.col(foundEigenValues) = tempEigenVectors.col(i);
                        foundEigenValues++;
                    }
                    else
                    {
                        std::cout<<"ZERO GARBAGE "<<tempEigenValues(i)<<std::endl;
                    }
                }
                omega+=1e-6;
                decompositionDone = true;
            }
        }
        catch(std::runtime_error e)
        {
            omega+=std::numeric_limits<double>::epsilon();
            std::cout<<"Increase omega"<<std::endl;
        }
    }
    vortBasisField = (curl*velBasisField);

    gravity = Eigen::VectorXd::Zero(decMesh.getNumFaces());
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0)))
                {
                    gravity(it->id) = -(it->normal.y>0?1:-1)*9.81;
                }
            }
        }
    }
    gravity = velBasisField.transpose()*gravity;

    std::vector<Vertex> vertices = mesh->getVertices();

    vorticityField = Eigen::VectorXd::Zero(decMesh.getNumEdges());
    //vorticityField.setRandom();
    //glm::uvec3 dims = decMesh.getDimensions();
    //vorticityField(decMesh.getNumEdges()/2) = 2e+64;
/*
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            glm::dvec3 p1=glm::dvec3(vertices[it->v1].pos);
            glm::dvec3 p2=glm::dvec3(vertices[it->v2].pos);
            if(std::abs(glm::dot(p1-p2,glm::dvec3(0.0,0.0,1.0)))>0.0)
            {
                //if((p1.z>=-0.3&&p1.z<=0.3&&p1.y==0.0&&p2.z>=-0.3&&p2.z<=0.3&&p2.y==0.0))
                //{
                //    vorticityField(it->id) = ((p1-p2).x>0?1.0:-1.0)*2;
                //}
                if((p1.x>=0.0&&p1.x<=0.0&&p1.y==0.0&&p2.x>=0.0&&p2.x<=0.0&&p2.y==0.0))
                {

                    vorticityField(it->id) = ((p1-p2).z>0?1.0:-1.0)*2;
                }
                //if((p1.x>=0.5&&p1.x<=0.7&&p1.y==0.0&&p2.x>=0.5&&p2.x<=0.7&&p2.y==0.0))
                //{
                //    vorticityField(it->id) = ((p1-p2).z>0?1.0:-1.0)*2;
                //}
                //else if(p1.x<=-0.5&&p1.x>=-0.7&&p1.y==0.0&&p2.x<=-0.5&&p2.x>=-0.7&&p2.y==0.0)
                //{
                //    vorticityField(it->id) = ((p1-p2).z>0?1.0:-1.0)*2;
                //}
            }

        }
    }
    setInitialVorticityField(vorticityField);*/


    velocityField = Eigen::VectorXd::Zero(decMesh.getNumFaces());
    glm::uvec3 dims = decMesh.getDimensions();
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                glm::dvec3 p1=glm::dvec3(vertices[it->v1].pos);
                glm::dvec3 p2=glm::dvec3(vertices[it->v2].pos);
                glm::dvec3 p3=glm::dvec3(vertices[it->v3].pos);
                glm::dvec3 p4=glm::dvec3(vertices[it->v4].pos);
                unsigned int zId = (it->id%(dims.x*dims.y+2*dims.x*dims.y-dims.x-dims.y));
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0))/*&&
                   (p1.y<=-0.7||p2.y<=-0.7||p3.y<=-0.7||p4.y<=-0.7)*/)
                        velocityField(it->id) = (it->normal.y>0?1:-1)*1.0;

            }
        }
    }
    setInitialVelocityField(velocityField);
}

void PSFSolver::buildAdvection()
{
    std::vector<Eigen::MatrixXd> wedges;
    wedges.resize(static_cast<unsigned int>(eigenValues.rows()));
    advection.resize(static_cast<unsigned int>(eigenValues.rows()));
    for(unsigned int i=0;i<velBasisField.cols();i++)
    {
        wedges[i] = Eigen::MatrixXd(decMesh.getNumEdges(),eigenValues.rows());
        wedges[i].setZero();
        advection[i] = Eigen::MatrixXd(decMesh.getNumEdges(),eigenValues.rows());
        advection[i].setZero();
    }
    std::vector<Vertex> vertices = mesh->getVertices();
    for(VoxelIterator fit=decMesh.getVoxelIteratorBegin();fit!=decMesh.getVoxelIteratorEnd();fit++)
    {
        if(fit->inside==GridState::INSIDE)
        {
            Face3D f1 = decMesh.getFace(fit->f1);
            Face3D f2 = decMesh.getFace(fit->f2);
            Face3D f3 = decMesh.getFace(fit->f3);
            Face3D f4 = decMesh.getFace(fit->f4);
            Face3D f5 = decMesh.getFace(fit->f5);
            Face3D f6 = decMesh.getFace(fit->f6);

            unsigned int ie1 = f1.e1;
            unsigned int ie2 = f1.e2;
            unsigned int ie3 = f1.e3;
            unsigned int ie4 = f1.e4;
            unsigned int ie5 = f2.e1;
            unsigned int ie6 = f2.e2;
            unsigned int ie7 = f2.e3;
            unsigned int ie8 = f2.e4;
            unsigned int ie9 = f5.e1;
            unsigned int ie10 = f5.e3;
            unsigned int ie11 = f6.e1;
            unsigned int ie12 = f6.e3;

            double s1 = decMesh.getFaceSignum(fit->f1,fit->v1,fit->v2,fit->v6,fit->v5);
            double s2 = decMesh.getFaceSignum(fit->f2,fit->v8,fit->v7,fit->v3,fit->v4);
            double s3 = decMesh.getFaceSignum(fit->f3,fit->v5,fit->v6,fit->v7,fit->v8);
            double s4 = decMesh.getFaceSignum(fit->f4,fit->v4,fit->v3,fit->v2,fit->v1);
            double s5 = decMesh.getFaceSignum(fit->f5,fit->v4,fit->v1,fit->v5,fit->v8);
            double s6 = decMesh.getFaceSignum(fit->f6,fit->v2,fit->v3,fit->v7,fit->v6);

            for(unsigned int i=0;i<nEigenFunctions;i++)
            {
                double vel1a = velBasisField(fit->f1,i);
                double vel2a = velBasisField(fit->f2,i);
                double vel3a = velBasisField(fit->f3,i);
                double vel4a = velBasisField(fit->f4,i);
                double vel5a = velBasisField(fit->f5,i);
                double vel6a = velBasisField(fit->f6,i);

                for(unsigned int j=0;j<nEigenFunctions;j++)
                {
                    double vel1b = velBasisField(fit->f1,j);
                    double vel2b = velBasisField(fit->f2,j);
                    double vel3b = velBasisField(fit->f3,j);
                    double vel4b = velBasisField(fit->f4,j);
                    double vel5b = velBasisField(fit->f5,j);
                    double vel6b = velBasisField(fit->f6,j);

                    wedges[i](ie1,j) += (0.25)*(vel1a*vel4b-vel1b*vel4a)*(glm::dot(glm::cross(f1.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie2,j) += (0.25)*(vel1a*vel6b-vel1b*vel6a)*(glm::dot(glm::cross(f1.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (0.25)*(vel1a*vel3b-vel1b*vel3a)*(glm::dot(glm::cross(f1.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (0.25)*(vel1a*vel5b-vel1b*vel5a)*(glm::dot(glm::cross(f1.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (0.25)*(vel2a*vel4b-vel2b*vel4a)*(glm::dot(glm::cross(f2.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (0.25)*(vel2a*vel5b-vel2b*vel5a)*(glm::dot(glm::cross(f2.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (0.25)*(vel2a*vel3b-vel2b*vel3a)*(glm::dot(glm::cross(f2.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (0.25)*(vel2a*vel6b-vel2b*vel6a)*(glm::dot(glm::cross(f2.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (0.25)*(vel5a*vel4b-vel5b*vel4a)*(glm::dot(glm::cross(f5.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (0.25)*(vel5a*vel3b-vel5b*vel3a)*(glm::dot(glm::cross(f5.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (0.25)*(vel6a*vel4b-vel6b*vel4a)*(glm::dot(glm::cross(f6.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (0.25)*(vel6a*vel3b-vel6b*vel3a)*(glm::dot(glm::cross(f6.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                }
            }
        }
    }

    for(unsigned int i=0;i<nEigenFunctions;i++)
    {
        for(unsigned int j=0;j<nEigenFunctions;j++)
        {
            advection[i].col(j) = eigenValues(i)*wedges[i].col(j);
        }
    }
    for(unsigned int i=0;i<eigenValues.rows();i++)
    {
        advection[i] = vortBasisField.transpose()*advection[i];
    }
}
