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
/*
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
    }*/


    velocityField = velBasisField*basisCoeff;

    glm::uvec3 dims = decMesh.getDimensions();
    std::vector<Vertex>& vertex = mesh->getVertices();
    #pragma omp parallel for
    for(std::vector<Particle>::iterator it=particles.begin();it<particles.end();it++)
    {
        if(simTime<it->lifeTime)
        {
            int yOfs = static_cast<int>((mesh->getAABB().min.y+mesh->getAABB().getExtent().y+it->position.y)/resolution);
            int xOfs = static_cast<int>((mesh->getAABB().min.x+mesh->getAABB().getExtent().x+it->position.x)/resolution);
            int zOfs = static_cast<int>((mesh->getAABB().min.z+mesh->getAABB().getExtent().z+it->position.z)/resolution);

            int zOfsMinus = static_cast<int>((mesh->getAABB().min.z+mesh->getAABB().getExtent().z+it->position.z-resolution/2)/resolution);
            int zOfsPlus = static_cast<int>((mesh->getAABB().min.z+mesh->getAABB().getExtent().z+it->position.z+resolution/2)/resolution);
            int yOfsMinus = static_cast<int>((mesh->getAABB().min.y+mesh->getAABB().getExtent().y+it->position.y-resolution/2)/resolution);
            int yOfsPlus = static_cast<int>((mesh->getAABB().min.y+mesh->getAABB().getExtent().y+it->position.y+resolution/2)/resolution);
            int xOfsMinus = static_cast<int>((mesh->getAABB().min.x+mesh->getAABB().getExtent().x+it->position.x-resolution/2)/resolution);
            int xOfsPlus = static_cast<int>((mesh->getAABB().min.x+mesh->getAABB().getExtent().x+it->position.x+resolution/2)/resolution);

            xOfs = glm::clamp(xOfs,0,int(dims.x-1));
            yOfs = glm::clamp(yOfs,0,int(dims.y-1));
            zOfs = glm::clamp(zOfs,0,int(dims.z-1));

            xOfsMinus = glm::clamp(xOfsMinus,0,int(dims.x-1));
            xOfsPlus = glm::clamp(xOfsPlus,0,int(dims.x-1));
            yOfsMinus = glm::clamp(yOfsMinus,0,int(dims.y-1));
            yOfsPlus = glm::clamp(yOfsPlus,0,int(dims.y-1));
            zOfsMinus = glm::clamp(zOfsMinus,0,int(dims.z-1));
            zOfsPlus = glm::clamp(zOfsPlus,0,int(dims.z-1));

            assert(zOfsMinus>=0);
            assert(zOfsPlus>=0);
            assert(yOfsMinus>=0);
            assert(yOfsPlus>=0);
            assert(xOfsMinus>=0);
            assert(xOfsPlus>=0);

            Voxel3D v1x = decMesh.getVoxel(decMesh.getVoxelIndex(xOfs,yOfsMinus,zOfsMinus));
            Voxel3D v2x = decMesh.getVoxel(decMesh.getVoxelIndex(xOfs,yOfsPlus,zOfsMinus));
            Voxel3D v3x = decMesh.getVoxel(decMesh.getVoxelIndex(xOfs,yOfsMinus,zOfsPlus));
            Voxel3D v4x = decMesh.getVoxel(decMesh.getVoxelIndex(xOfs,yOfsPlus,zOfsPlus));

            Voxel3D v1y = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsMinus,yOfs,zOfsMinus));
            Voxel3D v2y = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsPlus,yOfs,zOfsMinus));
            Voxel3D v3y = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsMinus,yOfs,zOfsPlus));
            Voxel3D v4y = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsPlus,yOfs,zOfsPlus));

            Voxel3D v1z = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsMinus,yOfsMinus,zOfs));
            Voxel3D v2z = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsPlus,yOfsMinus,zOfs));
            Voxel3D v3z = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsMinus,yOfsPlus,zOfs));
            Voxel3D v4z = decMesh.getVoxel(decMesh.getVoxelIndex(xOfsPlus,yOfsPlus,zOfs));

            Face3D xf = decMesh.getFace(v1x.f5);
            Face3D yf = decMesh.getFace(v1y.f3);
            Face3D zf = decMesh.getFace(v1z.f1);

            glm::vec3 cf1x = xf.center;
            glm::vec3 cf1y = yf.center;
            glm::vec3 cf1z = zf.center;

            double vel1x,vel2x,vel3x,vel4x,vel5x,vel6x,vel7x,vel8x;
            vel1x = decMesh.getFaceSignum(v1x.f5)*velocityField(decMesh.signedIdToIndex(v1x.f5));
            vel3x = -decMesh.getFaceSignum(v1x.f6)*velocityField(decMesh.signedIdToIndex(v1x.f6));
            vel2x = decMesh.getFaceSignum(v2x.f5)*velocityField(decMesh.signedIdToIndex(v2x.f5));
            vel4x = -decMesh.getFaceSignum(v2x.f6)*velocityField(decMesh.signedIdToIndex(v2x.f6));
            vel5x = decMesh.getFaceSignum(v3x.f5)*velocityField(decMesh.signedIdToIndex(v3x.f5));
            vel7x = -decMesh.getFaceSignum(v3x.f6)*velocityField(decMesh.signedIdToIndex(v3x.f6));
            vel6x = decMesh.getFaceSignum(v4x.f5)*velocityField(decMesh.signedIdToIndex(v4x.f5));
            vel8x = -decMesh.getFaceSignum(v4x.f6)*velocityField(decMesh.signedIdToIndex(v4x.f6));

            double vel1y,vel2y,vel3y,vel4y,vel5y,vel6y,vel7y,vel8y;
            vel1y = decMesh.getFaceSignum(v1y.f3)*velocityField(decMesh.signedIdToIndex(v1y.f3));
            vel2y = -decMesh.getFaceSignum(v1y.f4)*velocityField(decMesh.signedIdToIndex(v1y.f4));
            vel3y = decMesh.getFaceSignum(v2y.f3)*velocityField(decMesh.signedIdToIndex(v2y.f3));
            vel4y = -decMesh.getFaceSignum(v2y.f4)*velocityField(decMesh.signedIdToIndex(v2y.f4));
            vel5y = decMesh.getFaceSignum(v3y.f3)*velocityField(decMesh.signedIdToIndex(v3y.f3));
            vel6y = -decMesh.getFaceSignum(v3y.f4)*velocityField(decMesh.signedIdToIndex(v3y.f4));
            vel7y = decMesh.getFaceSignum(v4y.f3)*velocityField(decMesh.signedIdToIndex(v4y.f3));
            vel8y = -decMesh.getFaceSignum(v4y.f4)*velocityField(decMesh.signedIdToIndex(v4y.f4));

            double vel1z,vel2z,vel3z,vel4z,vel5z,vel6z,vel7z,vel8z;
            vel1z = decMesh.getFaceSignum(v1z.f1)*velocityField(decMesh.signedIdToIndex(v1z.f1));
            vel2z = decMesh.getFaceSignum(v3z.f1)*velocityField(decMesh.signedIdToIndex(v3z.f1));
            vel3z = decMesh.getFaceSignum(v2z.f1)*velocityField(decMesh.signedIdToIndex(v2z.f1));
            vel4z = decMesh.getFaceSignum(v4z.f1)*velocityField(decMesh.signedIdToIndex(v4z.f1));
            vel5z = -decMesh.getFaceSignum(v1z.f2)*velocityField(decMesh.signedIdToIndex(v1z.f2));
            vel6z = -decMesh.getFaceSignum(v3z.f2)*velocityField(decMesh.signedIdToIndex(v3z.f2));
            vel7z = -decMesh.getFaceSignum(v2z.f2)*velocityField(decMesh.signedIdToIndex(v2z.f2));
            vel8z = -decMesh.getFaceSignum(v4z.f2)*velocityField(decMesh.signedIdToIndex(v4z.f2));

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
    Eigen::SparseMatrix<double> mat2 = derivative1(decMesh,false)*derivative1(decMesh,true);
    //std::cout<<Eigen::MatrixXd(mat2)<<std::endl;
    //exit(1);

    Eigen::SparseMatrix<double> mat = -1.0*derivative1(decMesh,false)*hodge2(decMesh,true)*derivative1(decMesh,true)*hodge2(decMesh,false);
    //std::cout<<Eigen::MatrixXd(mat)<<std::endl;
    //exit(1);

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
    //mat.pruned();

    eigenValues.resize(nEigenFunctions);
    velBasisField.resize(decMesh.getNumFaces(),nEigenFunctions);
    bool decompositionDone=false;
    unsigned int foundEigenValues=0;
    double omega = 0.0;
    while(foundEigenValues<nEigenFunctions)
    {
        try
        {
            Spectra::SparseSymShiftSolve<double> op(mat);
            Spectra::SymEigsShiftSolver<double,Spectra::WHICH_LM,Spectra::SparseSymShiftSolve<double>> solver(&op,nEigenFunctions,2*nEigenFunctions,omega);
            solver.init();

            int nconv = solver.compute();//solver.compute(1000,1e-6,Spectra::WHICH_LM);
            if(solver.info()==Spectra::SUCCESSFUL)
            {

                Eigen::VectorXd tempEigenValues = solver.eigenvalues().real();
                Eigen::MatrixXd tempEigenVectors = solver.eigenvectors().real();
                bool onlyZero=true;
                for(unsigned int i=0;i<nconv && foundEigenValues<nEigenFunctions;i++)
                {
                    //if(std::abs(tempEigenValues(i))>0.0)
                    {
                        bool doubleEv = false;
                        for(unsigned int j=0;j<foundEigenValues;j++)
                        {
                            if(abs(velBasisField.col(j).dot(tempEigenVectors.col(i))-1.0)<std::numeric_limits<double>::epsilon())
                            {
                                std::cout<<"Already Inside"<<std::endl;
                                doubleEv = true;
                                break;
                            }
                        }
                        if(!doubleEv)
                        {
                        std::cout<<"ONE GARBAGE "<<tempEigenValues(i)<<std::endl;
                        //if(tempEigenValues(i)!=0.0)
                        {
                            omega = tempEigenValues(i);
                            onlyZero = false;
                        }
                        eigenValues(foundEigenValues) = tempEigenValues(i);
                        velBasisField.col(foundEigenValues) = tempEigenVectors.col(i);
                        foundEigenValues++;
                        }
                    }
                    //else
                    {
                        std::cout<<"ZERO GARBAGE "<<tempEigenValues(i)<<std::endl;
                    }
                }
                if(onlyZero)
                  omega+=1e-6;
                decompositionDone = true;
            }
            else
            {
                std::cout<<"Nothing Found"<<std::endl;
                omega+=1e-6;
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

                    vorticityField(it->id) = ((p1-p2).z>0?1.0:-1.0)*10;
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


    velocityField = Eigen::VectorXd::Zero(decMesh.getNumFaces());
    glm::uvec3 dims = decMesh.getDimensions();
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                AABB aabb = mesh->getAABB();
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0))&&
                   it->center.y<aabb.min.y+0.4)
                        velocityField(decMesh.getFaceIndex(*it)) = (it->normal.y>0?1:-1)*1.0;

            }
        }
    }
    setInitialVelocityField(velocityField);
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
    double scale = (decMesh.resolution/2)*(decMesh.resolution/2);

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

            unsigned int ie1 = decMesh.signedIdToIndex(f1.e1);
            unsigned int ie2 = decMesh.signedIdToIndex(f1.e2);
            unsigned int ie3 = decMesh.signedIdToIndex(f1.e3);
            unsigned int ie4 = decMesh.signedIdToIndex(f1.e4);
            unsigned int ie5 = decMesh.signedIdToIndex(f2.e1);
            unsigned int ie6 = decMesh.signedIdToIndex(f2.e2);
            unsigned int ie7 = decMesh.signedIdToIndex(f2.e3);
            unsigned int ie8 = decMesh.signedIdToIndex(f2.e4);
            unsigned int ie9 = decMesh.signedIdToIndex(f5.e1);
            unsigned int ie10 = decMesh.signedIdToIndex(f5.e3);
            unsigned int ie11 = decMesh.signedIdToIndex(f6.e1);
            unsigned int ie12 = decMesh.signedIdToIndex(f6.e3);

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

            Edge3D e1 = decMesh.getEdge(f1.e1);
            Edge3D e2 = decMesh.getEdge(f1.e2);
            Edge3D e3 = decMesh.getEdge(f1.e3);
            Edge3D e4 = decMesh.getEdge(f1.e4);
            Edge3D e5 = decMesh.getEdge(f2.e1);
            Edge3D e6 = decMesh.getEdge(f2.e2);
            Edge3D e7 = decMesh.getEdge(f2.e3);
            Edge3D e8 = decMesh.getEdge(f2.e4);
            Edge3D e9 = decMesh.getEdge(f5.e1);
            Edge3D e10 = decMesh.getEdge(f5.e3);
            Edge3D e11 = decMesh.getEdge(f6.e1);
            Edge3D e12 = decMesh.getEdge(f6.e3);

            Vertex e1v1 = vertices[e1.v1];
            Vertex e1v2 = vertices[e1.v2];
            Vertex e2v1 = vertices[e2.v1];
            Vertex e2v2 = vertices[e2.v2];
            Vertex e3v1 = vertices[e3.v1];
            Vertex e3v2 = vertices[e3.v2];
            Vertex e4v1 = vertices[e4.v1];
            Vertex e4v2 = vertices[e4.v2];

            Vertex e5v1 = vertices[e5.v1];
            Vertex e5v2 = vertices[e5.v2];
            Vertex e6v1 = vertices[e6.v1];
            Vertex e6v2 = vertices[e6.v2];
            Vertex e7v1 = vertices[e7.v1];
            Vertex e7v2 = vertices[e7.v2];
            Vertex e8v1 = vertices[e8.v1];
            Vertex e8v2 = vertices[e8.v2];

            glm::dvec3 evec1 = glm::dvec3(e1v2.pos-e1v1.pos);
            glm::dvec3 evec2 = glm::dvec3(e2v2.pos-e2v1.pos);
            glm::dvec3 evec3 = glm::dvec3(e3v2.pos-e3v1.pos);
            glm::dvec3 evec4 = glm::dvec3(e4v2.pos-e4v1.pos);

            glm::dvec3 evec5 = glm::dvec3(e5v2.pos-e5v1.pos);
            glm::dvec3 evec6 = glm::dvec3(e6v2.pos-e6v1.pos);
            glm::dvec3 evec7 = glm::dvec3(e7v2.pos-e7v1.pos);
            glm::dvec3 evec8 = glm::dvec3(e8v2.pos-e8v1.pos);

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
                    //assert(glm::dot(glm::cross(s1*f1.normal,-s4*f4.normal),evec1)>0.0);
                    //assert(glm::dot(glm::cross(-s1*f1.normal,s3*f3.normal),-evec3)>0.0);

                    //assert(glm::dot(glm::cross(s2*f2.normal,-s4*f4.normal),evec5)>0.0);
                    //assert(glm::dot(glm::cross(-s2*f2.normal,s3*f3.normal),evec7)>0.0);


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


                    wedges[i](ie1,j) += (scale)*(-vel1a*vel4b+vel1b*vel4a)*(glm::dot(glm::cross(-s1*f1.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie2,j) += (scale)*(-vel1a*vel6b+vel1b*vel6a)*(glm::dot(glm::cross(s1*f1.normal,-s6*f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (scale)*(-vel1a*vel3b+vel1b*vel3a)*(glm::dot(glm::cross(s1*f1.normal,-s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (scale)*(-vel1a*vel5b+vel1b*vel5a)*(glm::dot(glm::cross(-s1*f1.normal,s5*f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (scale)*(-vel2a*vel4b+vel2b*vel4a)*(glm::dot(glm::cross(s2*f2.normal,-s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (scale)*(-vel2a*vel5b+vel2b*vel5a)*(glm::dot(glm::cross(s2*f2.normal,-s5*f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (scale)*(-vel2a*vel3b+vel2b*vel3a)*(glm::dot(glm::cross(-s2*f2.normal,s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (scale)*(-vel2a*vel6b+vel2b*vel6a)*(glm::dot(glm::cross(-s2*f2.normal,s6*f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (scale)*(-vel5a*vel4b+vel5b*vel4a)*(glm::dot(glm::cross(-s5*f5.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (scale)*(-vel5a*vel3b+vel5b*vel3a)*(glm::dot(glm::cross(s5*f5.normal,-s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (scale)*(-vel6a*vel4b+vel6b*vel4a)*(glm::dot(glm::cross(-s6*f6.normal,s4*f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (scale)*(-vel6a*vel3b+vel6b*vel3a)*(glm::dot(glm::cross(s6*f6.normal,-s3*f3.normal),glm::dvec3(1.0,1.0,1.0)));
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
