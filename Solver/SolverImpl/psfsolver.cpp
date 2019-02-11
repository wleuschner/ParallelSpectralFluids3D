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

    velocityField = velBasisField*basisCoeff;

    /*#pragma omp parallel for
    for(std::vector<glm::dvec2>::iterator it=particles.begin();it<particles.end();it++)
    {
        unsigned int yOfs = static_cast<unsigned int>(it->y)/resolution;
        unsigned int xOfs = static_cast<unsigned int>(it->x)/resolution;

        unsigned int yOfsMinus = static_cast<unsigned int>(it->y-resolution/2)/resolution;
        unsigned int yOfsPlus = static_cast<unsigned int>(it->y+resolution/2)/resolution;
        unsigned int xOfsMinus = static_cast<unsigned int>(it->x-resolution/2)/resolution;
        unsigned int xOfsPlus = static_cast<unsigned int>(it->x+resolution/2)/resolution;

        Face2D f1x = decMesh.getFace(yOfsMinus*((mesh->getWidth()/resolution)+2)+xOfs);
        Face2D f2x = decMesh.getFace(yOfsPlus*((mesh->getWidth()/resolution)+2)+xOfs);
        Face2D f1y = decMesh.getFace(yOfs*((mesh->getWidth()/resolution)+2)+xOfsMinus);
        Face2D f2y = decMesh.getFace(yOfs*((mesh->getWidth()/resolution)+2)+xOfsPlus);


        glm::dvec2 cf1x = mesh->vertex[f1x.v1].pos+0.5*(mesh->vertex[f1x.v2].pos-mesh->vertex[f1x.v1].pos);
        glm::dvec2 cf2x = mesh->vertex[f2x.v1].pos+0.5*(mesh->vertex[f2x.v2].pos-mesh->vertex[f2x.v1].pos);
        glm::dvec2 cf1y = mesh->vertex[f1y.v4].pos+0.5*(mesh->vertex[f1y.v1].pos-mesh->vertex[f1y.v4].pos);
        glm::dvec2 cf2y = mesh->vertex[f2y.v4].pos+0.5*(mesh->vertex[f2y.v1].pos-mesh->vertex[f2y.v4].pos);

        double vel1x,vel2x,vel3x,vel4x;
        vel1x = decMesh.getEdgeSignum(f1x.e1,f1x.v1,f1x.v2)*velocityField(f1x.e1);
        vel3x = -decMesh.getEdgeSignum(f1x.e3,f1x.v3,f1x.v4)*velocityField(f1x.e3);
        vel2x = decMesh.getEdgeSignum(f2x.e1,f2x.v1,f2x.v2)*velocityField(f2x.e1);
        vel4x = -decMesh.getEdgeSignum(f2x.e3,f2x.v3,f2x.v4)*velocityField(f2x.e3);


        double vel1y,vel2y,vel3y,vel4y;
        vel1y = decMesh.getEdgeSignum(f1y.e4,f1y.v4,f1y.v1)*velocityField(f1y.e4);
        vel2y = -decMesh.getEdgeSignum(f1y.e2,f1y.v2,f1y.v3)*velocityField(f1y.e2);
        vel3y = decMesh.getEdgeSignum(f2y.e4,f2y.v4,f2y.v1)*velocityField(f2y.e4);
        vel4y = -decMesh.getEdgeSignum(f2y.e2,f2y.v2,f2y.v3)*velocityField(f2y.e2);


        glm::dvec2 particleNormalizedX;
        glm::dvec2 particleNormalizedY;
        glm::dvec2 vel;

        particleNormalizedX = (1.0/resolution)*((*it)-(cf1x));
        glm::dvec2 yVelXInterp = glm::mix(glm::dvec2(vel1x,vel3x),glm::dvec2(vel2x,vel4x),particleNormalizedX.y);
        vel.x = glm::mix(yVelXInterp.x,yVelXInterp.y,particleNormalizedX.x);

        particleNormalizedY = (1.0/resolution)*((*it)-(cf1y));
        glm::dvec2 yVelYInterp = glm::mix(glm::dvec2(vel1y,vel3y),glm::dvec2(vel2y,vel4y),particleNormalizedY.y);
        vel.y = glm::mix(yVelYInterp.x,yVelYInterp.y,particleNormalizedY.x);

        (*it) = (*it)+timeStep*vel;

    }*/
    vorticityField = curl*velocityField;
    maxRotation = vorticityField.cwiseAbs().maxCoeff();
    minRotation = vorticityField.minCoeff();
}

void PSFSolver::buildLaplace()
{
    Eigen::SparseMatrix<double> mat = derivative1(decMesh,false)*hodge2(decMesh,true)*derivative1(decMesh,true)*hodge2(decMesh,false);
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
            mat.prune([k](int i,int j,double v){return !(i==k||j==k);});
        }
    }

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
                            if((velBasisField.col(j)-tempEigenVectors.col(i)).norm()<std::numeric_limits<double>::epsilon())
                            {
                                std::cout<<"Already Inside"<<std::endl;
                                continue;
                            }
                        }
                        omega = tempEigenValues(i);
                        eigenValues(foundEigenValues) = tempEigenValues(i);
                        velBasisField.col(foundEigenValues) = tempEigenVectors.col(i);
                        foundEigenValues++;
                    }
                    else
                    {
                        std::cout<<"ZERO GARBAGE"<<std::endl;
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
    std::vector<Vertex> vertices = mesh->getVertices();

    /*
    vorticityField = Eigen::VectorXd::Zero(decMesh.getNumEdges());
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            glm::dvec2 point=mesh->vertex[it->id].pos;
            if((point.x==512&&
                (point.y==512||point.y==512)))
            {
                vorticityField(it->id) = 1000*2*3.141;
            }
        }
    }
    setInitialVorticityField(vorticityField);*/


    velocityField = Eigen::VectorXd::Zero(decMesh.getNumFaces());
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                velocityField(it->id) = (it->normal.y>0?1:-1)*0.1;
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
            double s2 = decMesh.getFaceSignum(fit->f2,fit->v3,fit->v4,fit->v8,fit->v7);
            double s3 = decMesh.getFaceSignum(fit->f3,fit->v7,fit->v8,fit->v5,fit->v6);
            double s4 = decMesh.getFaceSignum(fit->f4,fit->v4,fit->v3,fit->v2,fit->v1);
            double s5 = decMesh.getFaceSignum(fit->f5,fit->v4,fit->v1,fit->v5,fit->v8);
            double s6 = decMesh.getFaceSignum(fit->f6,fit->v2,fit->v3,fit->v7,fit->v6);

            for(unsigned int i=0;i<nEigenFunctions;i++)
            {
                double vel1a = s1*velBasisField(fit->f1,i);
                double vel2a = s2*velBasisField(fit->f2,i);
                double vel3a = s3*velBasisField(fit->f3,i);
                double vel4a = s4*velBasisField(fit->f4,i);
                double vel5a = s5*velBasisField(fit->f5,i);
                double vel6a = s6*velBasisField(fit->f6,i);

                for(unsigned int j=0;j<nEigenFunctions;j++)
                {
                    double vel1b = s1*velBasisField(fit->f1,j);
                    double vel2b = s2*velBasisField(fit->f2,j);
                    double vel3b = s3*velBasisField(fit->f3,j);
                    double vel4b = s4*velBasisField(fit->f4,j);
                    double vel5b = s5*velBasisField(fit->f5,j);
                    double vel6b = s6*velBasisField(fit->f6,j);

                    wedges[i](ie1,j) += (0.25)*(vel1a*vel4b-vel1b*vel4a)*(glm::dot(glm::cross(s1*f1.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie2,j) += (0.25)*(vel1a*vel6b-vel1b*vel6a)*(glm::dot(glm::cross(s1*f1.normal,s6*f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (0.25)*(vel1a*vel3b-vel1b*vel3a)*(glm::dot(glm::cross(s1*f1.normal,s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (0.25)*(vel1a*vel5b-vel1b*vel5a)*(glm::dot(glm::cross(s1*f1.normal,s5*f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (0.25)*(vel2a*vel4b-vel2b*vel4a)*(glm::dot(glm::cross(s2*f2.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (0.25)*(vel2a*vel5b-vel2b*vel5a)*(glm::dot(glm::cross(s2*f2.normal,s5*f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (0.25)*(vel2a*vel3b-vel2b*vel3a)*(glm::dot(glm::cross(s2*f2.normal,s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (0.25)*(vel2a*vel6b-vel2b*vel6a)*(glm::dot(glm::cross(s2*f2.normal,s6*f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (0.25)*(vel5a*vel4b-vel5b*vel4a)*(glm::dot(glm::cross(s5*f5.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (0.25)*(vel5a*vel3b-vel5b*vel3a)*(glm::dot(glm::cross(s5*f5.normal,s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (0.25)*(vel6a*vel4b-vel6b*vel4a)*(glm::dot(glm::cross(s6*f6.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (0.25)*(vel6a*vel3b-vel6b*vel3a)*(glm::dot(glm::cross(s6*f6.normal,s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
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
