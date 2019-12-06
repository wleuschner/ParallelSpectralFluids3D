#include"psfsolver.h"
#include <glm/gtc/matrix_transform.hpp>
#include<GL/glew.h>
#include<iostream>
#include<set>
#include<map>
#include"microprofile/microprofile.h"
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
    MICROPROFILE_SCOPEI("SolverCPU","integrate",0xFF000000);
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
        basisCoeff(k) *= std::exp(viscosity*eigenValues(k)*(timeStep));
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
            unsigned int yOfs = static_cast<unsigned int>((mesh->getAABB().min.y+mesh->getAABB().getExtent().y+it->position.y)/resolution)+1;
            unsigned int xOfs = static_cast<unsigned int>((mesh->getAABB().min.x+mesh->getAABB().getExtent().x+it->position.x)/resolution)+1;
            unsigned int zOfs = static_cast<unsigned int>((mesh->getAABB().min.z+mesh->getAABB().getExtent().z+it->position.z)/resolution)+1;

            unsigned int zOfsMinus = static_cast<unsigned int>((mesh->getAABB().min.z+mesh->getAABB().getExtent().z+it->position.z-resolution/2)/resolution)+1;
            unsigned int zOfsPlus = static_cast<unsigned int>((mesh->getAABB().min.z+mesh->getAABB().getExtent().z+it->position.z+resolution/2)/resolution)+1;
            unsigned int yOfsMinus = static_cast<unsigned int>((mesh->getAABB().min.y+mesh->getAABB().getExtent().y+it->position.y-resolution/2)/resolution)+1;
            unsigned int yOfsPlus = static_cast<unsigned int>((mesh->getAABB().min.y+mesh->getAABB().getExtent().y+it->position.y+resolution/2)/resolution)+1;
            unsigned int xOfsMinus = static_cast<unsigned int>((mesh->getAABB().min.x+mesh->getAABB().getExtent().x+it->position.x-resolution/2)/resolution)+1;
            unsigned int xOfsPlus = static_cast<unsigned int>((mesh->getAABB().min.x+mesh->getAABB().getExtent().x+it->position.x+resolution/2)/resolution)+1;

            assert(zOfsMinus>=0);
            assert(zOfsPlus>=0);
            assert(yOfsMinus>=0);
            assert(yOfsPlus>=0);
            assert(xOfsMinus>=0);
            assert(xOfsPlus>=0);

            Voxel3D v1x = decMesh.getVoxel((zOfsMinus*(dims.y*dims.x)+yOfsMinus*dims.x+xOfs)+1);
            Voxel3D v2x = decMesh.getVoxel((zOfsMinus*(dims.y*dims.x)+yOfsPlus*dims.x+xOfs)+1);
            Voxel3D v3x = decMesh.getVoxel((zOfsPlus*(dims.y*dims.x)+yOfsMinus*dims.x+xOfs)+1);
            Voxel3D v4x = decMesh.getVoxel((zOfsPlus*(dims.y*dims.x)+yOfsPlus*dims.x+xOfs)+1);

            Voxel3D v1y = decMesh.getVoxel((zOfsMinus*(dims.y*dims.x)+yOfs*dims.x+xOfsMinus)+1);
            Voxel3D v2y = decMesh.getVoxel((zOfsMinus*(dims.y*dims.x)+yOfs*dims.x+xOfsPlus)+1);
            Voxel3D v3y = decMesh.getVoxel((zOfsPlus*(dims.y*dims.x)+yOfs*dims.x+xOfsMinus)+1);
            Voxel3D v4y = decMesh.getVoxel((zOfsPlus*(dims.y*dims.x)+yOfs*dims.x+xOfsPlus)+1);

            Voxel3D v1z = decMesh.getVoxel((zOfs*(dims.y*dims.x)+yOfsMinus*dims.x+xOfsMinus)+1);
            Voxel3D v2z = decMesh.getVoxel((zOfs*(dims.y*dims.x)+yOfsMinus*dims.x+xOfsPlus)+1);
            Voxel3D v3z = decMesh.getVoxel((zOfs*(dims.y*dims.x)+yOfsPlus*dims.x+xOfsMinus)+1);
            Voxel3D v4z = decMesh.getVoxel((zOfs*(dims.y*dims.x)+yOfsPlus*dims.x+xOfsPlus)+1);

            Face3D xf = decMesh.getFace(v1x.f5);
            Face3D yf = decMesh.getFace(v1y.f3);
            Face3D zf = decMesh.getFace(v1z.f1);

            glm::vec3 cf1x = xf.center;
            glm::vec3 cf1y = yf.center;
            glm::vec3 cf1z = zf.center;
/*
            double vel1x,vel2x,vel3x,vel4x,vel5x,vel6x,vel7x,vel8x;
            vel1x = v1x.f5>0?1:-1*velocityField(labs(v1x.f5)-1);
            vel3x = v1x.f6>0?1:-1*velocityField(labs(v1x.f6)-1);
            vel2x = v2x.f5>0?1:-1*velocityField(labs(v2x.f5)-1);
            vel4x = v2x.f6>0?1:-1*velocityField(labs(v2x.f6)-1);
            vel5x = v3x.f5>0?1:-1*velocityField(labs(v3x.f5)-1);
            vel7x = v3x.f6>0?1:-1*velocityField(labs(v3x.f6)-1);
            vel6x = v4x.f5>0?1:-1*velocityField(labs(v4x.f5)-1);
            vel8x = v4x.f6>0?1:-1*velocityField(labs(v4x.f6)-1);

            double vel1y,vel2y,vel3y,vel4y,vel5y,vel6y,vel7y,vel8y;
            vel1y = v1y.f3>0?1:-1*velocityField(labs(v1y.f3)-1);
            vel2y = v1y.f4>0?1:-1*velocityField(labs(v1y.f4)-1);
            vel3y = v2y.f3>0?1:-1*velocityField(labs(v2y.f3)-1);
            vel4y = v2y.f4>0?1:-1*velocityField(labs(v2y.f4)-1);
            vel5y = v3y.f3>0?1:-1*velocityField(labs(v3y.f3)-1);
            vel6y = v3y.f4>0?1:-1*velocityField(labs(v3y.f4)-1);
            vel7y = v4y.f3>0?1:-1*velocityField(labs(v4y.f3)-1);
            vel8y = v4y.f4>0?1:-1*velocityField(labs(v4y.f4)-1);

            double vel1z,vel2z,vel3z,vel4z,vel5z,vel6z,vel7z,vel8z;
            vel1z = v1z.f1>0?1:-1*velocityField(labs(v1z.f1)-1);
            vel2z = v3z.f1>0?1:-1*velocityField(labs(v3z.f1)-1);
            vel3z = v2z.f1>0?1:-1*velocityField(labs(v2z.f1)-1);
            vel4z = v4z.f1>0?1:-1*velocityField(labs(v4z.f1)-1);
            vel5z = v1z.f2>0?1:-1*velocityField(labs(v1z.f2)-1);
            vel6z = v3z.f2>0?1:-1*velocityField(labs(v3z.f2)-1);
            vel7z = v2z.f2>0?1:-1*velocityField(labs(v2z.f2)-1);
            vel8z = v4z.f2>0?1:-1*velocityField(labs(v4z.f2)-1);*/

            double vel1x,vel2x,vel3x,vel4x,vel5x,vel6x,vel7x,vel8x;
            vel1x = velocityField(labs(v1x.f5)-1);
            vel3x = velocityField(labs(v1x.f6)-1);
            vel2x = velocityField(labs(v2x.f5)-1);
            vel4x = velocityField(labs(v2x.f6)-1);
            vel5x = velocityField(labs(v3x.f5)-1);
            vel7x = velocityField(labs(v3x.f6)-1);
            vel6x = velocityField(labs(v4x.f5)-1);
            vel8x = velocityField(labs(v4x.f6)-1);

            double vel1y,vel2y,vel3y,vel4y,vel5y,vel6y,vel7y,vel8y;
            vel1y = velocityField(labs(v1y.f3)-1);
            vel2y = velocityField(labs(v1y.f4)-1);
            vel3y = velocityField(labs(v2y.f3)-1);
            vel4y = velocityField(labs(v2y.f4)-1);
            vel5y = velocityField(labs(v3y.f3)-1);
            vel6y = velocityField(labs(v3y.f4)-1);
            vel7y = velocityField(labs(v4y.f3)-1);
            vel8y = velocityField(labs(v4y.f4)-1);

            double vel1z,vel2z,vel3z,vel4z,vel5z,vel6z,vel7z,vel8z;
            vel1z = velocityField(labs(v1z.f1)-1);
            vel2z = velocityField(labs(v3z.f1)-1);
            vel3z = velocityField(labs(v2z.f1)-1);
            vel4z = velocityField(labs(v4z.f1)-1);
            vel5z = velocityField(labs(v1z.f2)-1);
            vel6z = velocityField(labs(v3z.f2)-1);
            vel7z = velocityField(labs(v2z.f2)-1);
            vel8z = velocityField(labs(v4z.f2)-1);

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
    MICROPROFILE_SCOPEI("SolverCPU","buildLaplace",0xFF000000);
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
                    gravity(decMesh.getFaceIndex(*it)) = -(it->normal.y>0?1:-1)*0.00981;
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

/*
    velocityField = Eigen::VectorXd::Zero(decMesh.getNumFaces());
    glm::uvec3 dims = decMesh.getDimensions();
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                //glm::dvec3 p1=glm::dvec3(vertices[it->v1].pos);
                //glm::dvec3 p2=glm::dvec3(vertices[it->v2].pos);
                //glm::dvec3 p3=glm::dvec3(vertices[it->v3].pos);
                //glm::dvec3 p4=glm::dvec3(vertices[it->v4].pos);
                AABB aabb = mesh->getAABB();
                unsigned int zId = (decMesh.getFaceIndex(*it)%(dims.x*dims.y+2*dims.x*dims.y-dims.x-dims.y));
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0))&&
                   it->center.y<aabb.min.y+0.2)
                        velocityField(decMesh.getFaceIndex(*it)) = (it->normal.y>0?1:-1)*0.1;

            }
        }
    }
    setInitialVelocityField(velocityField);*/
}

void PSFSolver::buildAdvection()
{
    MICROPROFILE_SCOPEI("SolverCPU","buildAdvection",MP_RED);
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

            unsigned int ie1 = labs(f1.e1)-1;
            unsigned int ie2 = labs(f1.e2)-1;
            unsigned int ie3 = labs(f1.e3)-1;
            unsigned int ie4 = labs(f1.e4)-1;
            unsigned int ie5 = labs(f2.e1)-1;
            unsigned int ie6 = labs(f2.e2)-1;
            unsigned int ie7 = labs(f2.e3)-1;
            unsigned int ie8 = labs(f2.e4)-1;
            unsigned int ie9 = labs(f5.e1)-1;
            unsigned int ie10 = labs(f5.e3)-1;
            unsigned int ie11 = labs(f6.e1)-1;
            unsigned int ie12 = labs(f6.e3)-1;

            double se1 = f1.e1>0?1:-1;
            double se2 = f1.e2>0?1:-1;
            double se3 = f1.e3>0?1:-1;
            double se4 = f1.e4>0?1:-1;
            double se5 = f2.e1>0?1:-1;
            double se6 = f2.e2>0?1:-1;
            double se7 = f2.e3>0?1:-1;
            double se8 = f2.e4>0?1:-1;
            double se9 = f5.e1>0?1:-1;
            double se10 = f5.e3>0?1:-1;
            double se11 = f6.e1>0?1:-1;
            double se12 = f6.e3>0?1:-1;
/*
            Edge3D e1 = decMesh.getEdge(f1.e1);
            Edge3D e2 = decMesh.getEdge(f1.e2);
            Edge3D e3 = decMesh.getEdge(f1.e3);
            Edge3D e4 = decMesh.getEdge(f1.e4);

            Vertex e1v1 = vertices[e1.v1];
            Vertex e1v2 = vertices[e1.v2];
            Vertex e2v1 = vertices[e2.v1];
            Vertex e2v2 = vertices[e2.v2];
            Vertex e3v1 = vertices[e3.v1];
            Vertex e3v2 = vertices[e3.v2];
            Vertex e4v1 = vertices[e4.v1];
            Vertex e4v2 = vertices[e4.v2];

            glm::dvec3 evec1 = glm::dvec3(e1v2.pos-e1v1.pos);
            glm::dvec3 evec2 = glm::dvec3(e2v2.pos-e2v1.pos);
            glm::dvec3 evec3 = glm::dvec3(e3v2.pos-e3v1.pos);
            glm::dvec3 evec4 = glm::dvec3(e4v2.pos-e4v1.pos);*/

            double s1 = fit->f1>0?1:-1;
            double s2 = fit->f2>0?1:-1;
            double s3 = fit->f3>0?1:-1;
            double s4 = fit->f4>0?1:-1;
            double s5 = fit->f5>0?1:-1;
            double s6 = fit->f6>0?1:-1;

            for(unsigned int i=0;i<nEigenFunctions;i++)
            {
                double vel1a = s1*velBasisField(labs(fit->f1)-1,i);
                double vel2a = s2*velBasisField(labs(fit->f2)-1,i);
                double vel3a = s3*velBasisField(labs(fit->f3)-1,i);
                double vel4a = s4*velBasisField(labs(fit->f4)-1,i);
                double vel5a = s5*velBasisField(labs(fit->f5)-1,i);
                double vel6a = s6*velBasisField(labs(fit->f6)-1,i);

                for(unsigned int j=0;j<nEigenFunctions;j++)
                {
                    double vel1b = s1*velBasisField(labs(fit->f1)-1,j);
                    double vel2b = s2*velBasisField(labs(fit->f2)-1,j);
                    double vel3b = s3*velBasisField(labs(fit->f3)-1,j);
                    double vel4b = s4*velBasisField(labs(fit->f4)-1,j);
                    double vel5b = s5*velBasisField(labs(fit->f5)-1,j);
                    double vel6b = s6*velBasisField(labs(fit->f6)-1,j);

                    //assert(glm::dot(evec1,glm::cross(f1.normal,f4.normal))>0.0);
                    //assert(glm::dot(evec1,se2*glm::cross(s1*f1.normal,-s6*f6.normal))>0.0);
                    //assert(glm::dot(evec1,se3*glm::cross(s1*f1.normal,-s3*f3.normal))>0.0);
                    //assert(glm::dot(evec1,se4*glm::cross(-s1*f1.normal,s5*f5.normal))>0.0);


                    wedges[i](ie1,j) += (0.25)*(-vel1a*vel4b+vel1b*vel4a)*(glm::dot(glm::cross(-s1*f1.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie2,j) += (0.25)*(-vel1a*vel6b+vel1b*vel6a)*(glm::dot(glm::cross(s1*f1.normal,-s6*f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (0.25)*(-vel1a*vel3b+vel1b*vel3a)*(glm::dot(glm::cross(s1*f1.normal,-s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (0.25)*(-vel1a*vel5b+vel1b*vel5a)*(glm::dot(glm::cross(-s1*f1.normal,s5*f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (0.25)*(-vel2a*vel4b+vel2b*vel4a)*(glm::dot(glm::cross(s2*f2.normal,-s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (0.25)*(-vel2a*vel5b+vel2b*vel5a)*(glm::dot(glm::cross(s2*f2.normal,-s5*f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (0.25)*(-vel2a*vel3b+vel2b*vel3a)*(glm::dot(glm::cross(-s2*f2.normal,s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (0.25)*(-vel2a*vel6b+vel2b*vel6a)*(glm::dot(glm::cross(-s2*f2.normal,s6*f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (0.25)*(-vel5a*vel4b+vel5b*vel4a)*(glm::dot(glm::cross(-s5*f5.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (0.25)*(-vel5a*vel3b+vel5b*vel3a)*(glm::dot(glm::cross(s5*f5.normal,-s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (0.25)*(-vel6a*vel4b+vel6b*vel4a)*(glm::dot(glm::cross(-s6*f6.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (0.25)*(-vel6a*vel3b+vel6b*vel3a)*(glm::dot(glm::cross(s6*f6.normal,-s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                }
            }
        }
    }

    #pragma omp parallel for
    for(unsigned int i=0;i<nEigenFunctions;i++)
    {
        for(unsigned int j=0;j<nEigenFunctions;j++)
        {
            advection[i].col(j) = eigenValues(i)*wedges[i].col(j);
        }
    }

    #pragma omp parallel for
    for(unsigned int i=0;i<eigenValues.rows();i++)
    {
        advection[i] = vortBasisField.transpose()*advection[i];
    }
}
