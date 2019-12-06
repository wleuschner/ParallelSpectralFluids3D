#include"psfsolvergpu.h"
#include <GL/glew.h>
#include <viennacl/scalar.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/ocl/backend.hpp>
#include <fstream>
#include <string>
#include <CL/cl_gl.h>
#include "../../Spectra/MatOp/SparseSymShiftSolve.h"
#include "../../Spectra/SymEigsShiftSolver.h"
#include "../../DEC/dec.h"

PSFSolverGPU::PSFSolverGPU() : AbstractSolver()
{
}

void PSFSolverGPU::integrate()
{
    viennacl::scalar<double> e1 = 0.0;
    viennacl::scalar<double> e2 = 0.0;

    e1 = viennacl::linalg::inner_prod(vcl_basisCoeff,vcl_basisCoeff);

    viennacl::vector<double> vel = viennacl::zero_vector<double>(nEigenFunctions);

    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        vel(k) += viennacl::linalg::inner_prod(vcl_basisCoeff,viennacl::linalg::prod(vcl_advection[k],vcl_basisCoeff));
    }
    vcl_basisCoeff += timeStep*vel;

    e2 = viennacl::linalg::inner_prod(vcl_basisCoeff,vcl_basisCoeff);

    basisCoeff *= std::sqrt(e1/e2);

    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        vcl_basisCoeff(k) *= std::exp(-viscosity*eigenValues(k)*(timeStep));
    }
    if(gravityActive)
    {
        vcl_basisCoeff += timeStep*vcl_gravity;
    }
    vcl_velocityField = viennacl::linalg::prod(vcl_velBasisField,vcl_basisCoeff);
}

void PSFSolverGPU::buildLaplace()
{
    std::ifstream f("Resources/Simulation/simulation.cl");
    std::string source((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());

//    cl_mem particlesBuffer = clCreateFromGLBuffer(viennacl::ocl::current_context(),CL_MEM_READ_WRITE,particleVerts->id);
    viennacl::ocl::program& vcl_prog_psf = viennacl::ocl::current_context().add_program(source.data(),"advection");
    viennacl::ocl::kernel& advection_kernel = vcl_prog_psf.get_kernel("advection");

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
            mat.prune([k](int i,int j,double v){return !(i==k||j==k);});
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
                    gravity(it->id) = -(it->normal.y>0?1:-1)*0.00981;
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
    setInitialVorticityField(vorticityField);
    */

    velocityField = Eigen::VectorXd::Zero(decMesh.getNumFaces());
    glm::uvec3 dims = decMesh.getDimensions();
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                AABB aabb = mesh->getAABB();
                unsigned int zId = (it->id%(dims.x*dims.y+2*dims.x*dims.y-dims.x-dims.y));
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0))&&
                   it->center.y<=aabb.min.y+0.1)
                        velocityField(it->id) = (it->normal.y>0?1:-1)*0.00001;

            }
        }
    }
    setInitialVelocityField(velocityField);
    vcl_vortBasisField.resize(vortBasisField.rows(),vortBasisField.cols());
    vcl_velBasisField.resize(velBasisField.rows(),velBasisField.cols());
    vcl_gravity.resize(gravity.size());
    vcl_vorticityField.resize(vorticityField.rows(),vorticityField.cols());
    vcl_velocityField.resize(velBasisField.rows(),velBasisField.cols());
    vcl_basisCoeff.resize(basisCoeff.size());

    viennacl::copy(vortBasisField,vcl_vortBasisField);
    viennacl::copy(velBasisField,vcl_velBasisField);
    viennacl::copy(gravity,vcl_gravity);
    viennacl::copy(vorticityField,vcl_vorticityField);
    viennacl::copy(velocityField,vcl_velocityField);
    viennacl::copy(basisCoeff,vcl_basisCoeff);
}

void PSFSolverGPU::buildAdvection()
{
    std::vector<Eigen::MatrixXd> wedges;
    wedges.resize(static_cast<unsigned int>(eigenValues.rows()));
    advection.resize(static_cast<unsigned int>(eigenValues.rows()));
    vcl_advection.resize(static_cast<unsigned int>(eigenValues.rows()));
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
        vcl_advection[i].resize(advection[i].rows(),advection[i].cols());
        viennacl::copy(advection[i],vcl_advection[i]);
    }
}
