#include"psfsolver.h"
#include <glm/gtc/matrix_transform.hpp>
#include<GL/glew.h>
#include<iostream>
#include<set>
#include<map>
#include<chrono>
#include"microprofile/microprofile.h"
#include "../../Spectra/MatOp/SparseGenRealShiftSolve.h"
#include "../../Spectra/GenEigsRealShiftSolver.h"
#include"../../Graphics/TextureArray/TextureArray.h"
#include"../../Graphics/FrameBufferObject/FrameBufferObject.h"
#include"../../DEC/dec.h"

PSFSolver::PSFSolver() : AbstractSolver()
{
    gpu=false;
    simTime = 0.0f;
}

void PSFSolver::integrate()
{
    if(benchmarkFrameNo>maxFramesBenchmark)
    {
        startBenchmark(false);
    }
    if(isBenchmark)
    {
        benchmarkFrameNo++;
    }
    {
        std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();

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
        velocityField = (velBasisField*basisCoeff).cast<float>();
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        if(isBenchmark)
        {
            double t = std::chrono::duration<double,std::milli>(end-begin).count();
            benchmarkSums[0] +=t;
            benchmark_file<<t<<",";
        }
    }

    {
        std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
        glm::uvec3 dims = decMesh.getDimensions();
        std::vector<Vertex>& vertex = mesh->getVertices();
        std::vector<Particle>& parts = particles->getParticles();
        #pragma omp parallel for
        for(std::vector<Particle>::iterator it=parts.begin();it<parts.end();it++)
        {
            if(it->position.w>0.0f)
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

                unsigned int v1xIdx = decMesh.getVoxelIndex(xOfs,yOfsMinus,zOfsMinus);
                unsigned int v2xIdx = decMesh.getVoxelIndex(xOfs,yOfsPlus,zOfsMinus);
                unsigned int v3xIdx = decMesh.getVoxelIndex(xOfs,yOfsMinus,zOfsPlus);
                unsigned int v4xIdx = decMesh.getVoxelIndex(xOfs,yOfsPlus,zOfsPlus);

                unsigned int v1yIdx = decMesh.getVoxelIndex(xOfsMinus,yOfs,zOfsMinus);
                unsigned int v2yIdx = decMesh.getVoxelIndex(xOfsPlus,yOfs,zOfsMinus);
                unsigned int v3yIdx = decMesh.getVoxelIndex(xOfsMinus,yOfs,zOfsPlus);
                unsigned int v4yIdx = decMesh.getVoxelIndex(xOfsPlus,yOfs,zOfsPlus);

                unsigned int v1zIdx = decMesh.getVoxelIndex(xOfsMinus,yOfsMinus,zOfs);
                unsigned int v2zIdx = decMesh.getVoxelIndex(xOfsPlus,yOfsMinus,zOfs);
                unsigned int v3zIdx = decMesh.getVoxelIndex(xOfsMinus,yOfsPlus,zOfs);
                unsigned int v4zIdx = decMesh.getVoxelIndex(xOfsPlus,yOfsPlus,zOfs);

                glm::vec3 cf1x = glm::vec3(mesh->getAABB().min)+(float)(resolution)*glm::vec3(xOfs,yOfsMinus,zOfsMinus)+0.5f*glm::vec3(0.0f,resolution,resolution);

                glm::vec3 cf1y = glm::vec3(mesh->getAABB().min)+(float)(resolution)*glm::vec3(xOfsMinus,yOfs,zOfsMinus)+0.5f*glm::vec3(resolution,0.0f,resolution);

                glm::vec3 cf1z = glm::vec3(mesh->getAABB().min)+(float)(resolution)*glm::vec3(xOfsMinus,yOfsMinus,zOfs)+0.5f*glm::vec3(resolution,resolution,0.0f);

                unsigned int v1xf5 = decMesh.getXFaceIndex(xOfs,yOfsMinus,zOfsMinus);
                unsigned int v1xf6 = decMesh.getXFaceIndex(xOfs+1,yOfsMinus,zOfsMinus);
                unsigned int v2xf5 = decMesh.getXFaceIndex(xOfs,yOfsPlus,zOfsMinus);
                unsigned int v2xf6 = decMesh.getXFaceIndex(xOfs+1,yOfsPlus,zOfsMinus);
                unsigned int v3xf5 = decMesh.getXFaceIndex(xOfs,yOfsMinus,zOfsPlus);
                unsigned int v3xf6 = decMesh.getXFaceIndex(xOfs+1,yOfsMinus,zOfsPlus);
                unsigned int v4xf5 = decMesh.getXFaceIndex(xOfs,yOfsPlus,zOfsPlus);
                unsigned int v4xf6 = decMesh.getXFaceIndex(xOfs+1,yOfsPlus,zOfsPlus);


                unsigned int v1yf3 = decMesh.getYFaceIndex(xOfsMinus,yOfs,zOfsMinus);
                unsigned int v1yf4 = decMesh.getYFaceIndex(xOfsMinus,yOfs+1,zOfsMinus);
                unsigned int v2yf3 = decMesh.getYFaceIndex(xOfsPlus,yOfs,zOfsMinus);
                unsigned int v2yf4 = decMesh.getYFaceIndex(xOfsPlus,yOfs+1,zOfsMinus);
                unsigned int v3yf3 = decMesh.getYFaceIndex(xOfsMinus,yOfs,zOfsPlus);
                unsigned int v3yf4 = decMesh.getYFaceIndex(xOfsMinus,yOfs+1,zOfsPlus);
                unsigned int v4yf3 = decMesh.getYFaceIndex(xOfsPlus,yOfs,zOfsPlus);
                unsigned int v4yf4 = decMesh.getYFaceIndex(xOfsPlus,yOfs+1,zOfsPlus);


                unsigned int v1zf1 = decMesh.getZFaceIndex(xOfsMinus,yOfsMinus,zOfs);
                unsigned int v1zf2 = decMesh.getZFaceIndex(xOfsMinus,yOfsMinus,zOfs+1);
                unsigned int v2zf1 = decMesh.getZFaceIndex(xOfsPlus,yOfsMinus,zOfs);
                unsigned int v2zf2 = decMesh.getZFaceIndex(xOfsPlus,yOfsMinus,zOfs+1);
                unsigned int v3zf1 = decMesh.getZFaceIndex(xOfsMinus,yOfsPlus,zOfs);
                unsigned int v3zf2 = decMesh.getZFaceIndex(xOfsMinus,yOfsPlus,zOfs+1);
                unsigned int v4zf1 = decMesh.getZFaceIndex(xOfsPlus,yOfsPlus,zOfs);
                unsigned int v4zf2 = decMesh.getZFaceIndex(xOfsPlus,yOfsPlus,zOfs+1);


                float vel1x,vel2x,vel3x,vel4x,vel5x,vel6x,vel7x,vel8x;
                vel1x = -decMesh.getFaceSignum(v1xIdx,4)*velocityField((v1xf5));
                vel3x = decMesh.getFaceSignum(v1xIdx,5)*velocityField((v1xf6));
                vel2x = -decMesh.getFaceSignum(v2xIdx,4)*velocityField((v2xf5));
                vel4x = decMesh.getFaceSignum(v2xIdx,5)*velocityField((v2xf6));
                vel5x = -decMesh.getFaceSignum(v3xIdx,4)*velocityField((v3xf5));
                vel7x = decMesh.getFaceSignum(v3xIdx,5)*velocityField((v3xf6));
                vel6x = -decMesh.getFaceSignum(v4xIdx,4)*velocityField((v4xf5));
                vel8x = decMesh.getFaceSignum(v4xIdx,5)*velocityField((v4xf6));

                float vel1y,vel2y,vel3y,vel4y,vel5y,vel6y,vel7y,vel8y;
                vel1y = -decMesh.getFaceSignum(v1yIdx,2)*velocityField((v1yf3));
                vel2y = decMesh.getFaceSignum(v1yIdx,3)*velocityField((v1yf4));
                vel3y = -decMesh.getFaceSignum(v2yIdx,2)*velocityField((v2yf3));
                vel4y = decMesh.getFaceSignum(v2yIdx,3)*velocityField((v2yf4));
                vel5y = -decMesh.getFaceSignum(v3yIdx,2)*velocityField((v3yf3));
                vel6y = decMesh.getFaceSignum(v3yIdx,3)*velocityField((v3yf4));
                vel7y = -decMesh.getFaceSignum(v4yIdx,2)*velocityField((v4yf3));
                vel8y = decMesh.getFaceSignum(v4yIdx,3)*velocityField((v4yf4));

                float vel1z,vel2z,vel3z,vel4z,vel5z,vel6z,vel7z,vel8z;
                vel1z = decMesh.getFaceSignum(v1zIdx,0)*velocityField((v1zf1));
                vel2z = decMesh.getFaceSignum(v3zIdx,0)*velocityField((v3zf1));
                vel3z = decMesh.getFaceSignum(v2zIdx,0)*velocityField((v2zf1));
                vel4z = decMesh.getFaceSignum(v4zIdx,0)*velocityField((v4zf1));
                vel5z = -decMesh.getFaceSignum(v1zIdx,1)*velocityField((v1zf2));
                vel6z = -decMesh.getFaceSignum(v3zIdx,1)*velocityField((v3zf2));
                vel7z = -decMesh.getFaceSignum(v2zIdx,1)*velocityField((v2zf2));
                vel8z = -decMesh.getFaceSignum(v4zIdx,1)*velocityField((v4zf2));

                glm::vec3 particleNormalizedX;
                glm::vec3 particleNormalizedY;
                glm::vec3 particleNormalizedZ;
                glm::vec3 vel;

                particleNormalizedX = float(1.0/resolution)*((glm::vec3(it->position))-glm::vec3(cf1x));


                glm::vec4 xVelXInterp = glm::mix(glm::vec4(vel1x,vel3x,vel5x,vel7x),glm::vec4(vel2x,vel4x,vel6x,vel8x),particleNormalizedX.y);
                glm::vec2 xVelYInterp = glm::mix(glm::vec2(xVelXInterp.x,xVelXInterp.z),glm::vec2(xVelXInterp.y,xVelXInterp.w),particleNormalizedX.x);
                vel.x = glm::mix(xVelYInterp.x,xVelYInterp.y,particleNormalizedX.z);

                particleNormalizedY = float(1.0/resolution)*(glm::vec3(it->position)-glm::vec3(cf1y));

                glm::vec4 yVelXInterp = glm::mix(glm::vec4(vel1y,vel3y,vel5y,vel7y),glm::vec4(vel2y,vel4y,vel6y,vel8y),particleNormalizedY.y);
                glm::vec2 yVelYInterp = glm::mix(glm::vec2(yVelXInterp.x,yVelXInterp.z),glm::vec2(yVelXInterp.y,yVelXInterp.w),particleNormalizedY.x);
                vel.y = glm::mix(yVelYInterp.x,yVelYInterp.y,particleNormalizedY.z);

                particleNormalizedZ = float(1.0/resolution)*(glm::vec3(it->position)-glm::vec3(cf1z));

                glm::vec4 zVelXInterp = glm::mix(glm::vec4(vel1z,vel3z,vel5z,vel7z),glm::vec4(vel2z,vel4z,vel6z,vel8z),particleNormalizedZ.y);
                glm::vec2 zVelYInterp = glm::mix(glm::vec2(zVelXInterp.x,zVelXInterp.z),glm::vec2(zVelXInterp.y,zVelXInterp.w),particleNormalizedZ.x);
                vel.z = glm::mix(zVelYInterp.x,zVelYInterp.y,particleNormalizedZ.z);

                //it->position = (it->position)+timeStep*vel;
                it->position = glm::vec4(glm::clamp(glm::vec3(it->position)+float(timeStep)*vel,glm::vec3(mesh->getAABB().min),glm::vec3(mesh->getAABB().max)),it->position.w-1.0f);
            }
            else
            {
                it->position.x = ((rand()%1024)/1024.0-0.5)*0.25;
                it->position.y = mesh->getAABB().min.y+0.2f;
                it->position.z = -glm::dvec3(mesh->getAABB().getCenter()).z+((rand()%1024)/1024.0-0.5)*0.25;
                it->position.w = lifeTime*((rand()%1024)/1024.0);
            }
        }
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        if(isBenchmark)
        {
            double t = std::chrono::duration<double,std::milli>(end-begin).count();
            benchmarkSums[1] +=t;
            benchmark_file<<t<<std::endl;
        }
    }
    simTime+=1.0;
}

void PSFSolver::buildLaplace()
{
    MICROPROFILE_SCOPEI("SolverCPU","buildLaplace",0xFF000000);
    //Eigen::SparseMatrix<double> mat2 = derivative1(decMesh,false)*derivative1(decMesh,true);
    //std::cout<<Eigen::MatrixXd(derivative1(decMesh,true))<<std::endl;

    Eigen::SparseMatrix<double> mat3 = hodge2(decMesh,true);
    Eigen::SparseMatrix<double> mat4 = hodge2(decMesh,false);
    //std::cout<<Eigen::MatrixXd(mat3)<<std::endl;
    //std::cout<<Eigen::MatrixXd(mat4)<<std::endl;

    Eigen::SparseMatrix<double> mat = 1.0*derivative1(decMesh,false)*hodge2(decMesh,true)*derivative1(decMesh,true)*hodge2(decMesh,false);
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
            mat.prune([k](int i,int j,double v){return !(i==k);});
        }
    }
    //mat.pruned();

    eigenValues.resize(nEigenFunctions);
    velBasisField.resize(decMesh.getNumFaces(),nEigenFunctions);
    bool decompositionDone=false;
    unsigned int foundEigenValues=0;
    double omega = 0.1;
    while(foundEigenValues<nEigenFunctions)
    {
        try
        {
            Spectra::SparseGenRealShiftSolve<double> op(mat);
            Spectra::GenEigsRealShiftSolver<double,Spectra::WHICH_LM,Spectra::SparseGenRealShiftSolve<double>> solver(&op,nEigenFunctions,2*nEigenFunctions,omega);
            solver.init();

            int nconv = solver.compute(1000,1e-10,Spectra::WHICH_SM);
            if(solver.info()==Spectra::SUCCESSFUL)
            {
                Eigen::VectorXcd tempEigenValues = solver.eigenvalues();
                Eigen::MatrixXcd tempEigenVectors = solver.eigenvectors();
                bool onlyZero=true;
                for(unsigned int i=0;i<nconv && foundEigenValues<nEigenFunctions;i++)
                {
                    if(tempEigenValues(i).real()<0) continue;
                    if((std::abs(1.0/tempEigenValues(i).real())<1e-10 || std::abs(1.0/tempEigenValues(i).real())>1e+10)) continue;

                    Eigen::VectorXd backProjection = mat*tempEigenVectors.col(i).real();
                    if(backProjection.isApprox(tempEigenValues(i).real()*tempEigenVectors.col(i).real(),1))
                    //if(std::abs(1.0/tempEigenValues(i))>1e-10 && std::abs(1.0/tempEigenValues(i))<1e+10)
                    {
                        bool doubleEv = false;

                        /*for(unsigned int j=0;j<foundEigenValues;j++)
                        {
                            if(fabs(eigenValues(j)-tempEigenValues(i))<=0.1)
                            //if(velBasisField.col(j).isApprox(tempEigenVectors.col(i)))
                            {
                                std::cout<<"Already Inside "<<tempEigenValues(i)<<std::endl;
                                doubleEv = true;
                                break;
                            }
                        }*/
                        if(!doubleEv)
                        {
                        std::cout<<"ONE GARBAGE "<<tempEigenValues(i)<<std::endl;
                        //if(tempEigenValues(i)!=0.0)
                        {
                            omega = tempEigenValues(i).real();
                            onlyZero = false;
                        }
                        eigenValues(foundEigenValues) = tempEigenValues(i).real();
                        velBasisField.col(foundEigenValues) = tempEigenVectors.col(i).real();
                        foundEigenValues++;
                        }
                    }
                    else
                    {
                        std::cout<<"ZERO GARBAGE "<<tempEigenValues(i).real()<<std::endl;
                    }
                }
                if(onlyZero)
                  omega+=0.1;
                decompositionDone = true;
            }
            else
            {
                std::cout<<"Nothing Found"<<std::endl;
                omega+=0.1;
            }
        }
        catch(std::runtime_error e)
        {
            omega+=0.1;
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
                    gravity(decMesh.getFaceIndex(*it)) = -(it->normal.y>0?-1:1)*0.00981;
                }
            }
        }
    }
    gravity = velBasisField.transpose()*gravity;

    std::vector<Vertex> vertices = mesh->getVertices();

    vorticityField = Eigen::VectorXf::Zero(decMesh.getNumEdges());
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

                    vorticityField(it->id) = ((p1-p2).z>0?1.0:-1.0)*1;
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
    //setInitialVorticityField(vorticityField);

    velocityField = Eigen::VectorXf::Zero(decMesh.getNumFaces());
    glm::uvec3 dims = decMesh.getDimensions();
    for(FaceIterator it=decMesh.getFaceIteratorBegin();it!=decMesh.getFaceIteratorEnd();it++)
    {
        if(it->inside==GridState::INSIDE)
        {
            if(std::abs(glm::dot(glm::dvec3(0.0,1.0,0.0),it->normal))>std::numeric_limits<double>::epsilon())
            {
                AABB aabb = mesh->getAABB();
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0)) &&
                   it->center.y<aabb.min.y+0.8)
                        velocityField(decMesh.getFaceIndex(*it)) = (it->normal.y>0?-1:1)*1.0;

            }
        }
    }
    setInitialVelocityField(velocityField);

    for(unsigned int i=0;i<particles->getNumParticles();i++)
    {
        glm::dvec3 pos = glm::dvec3(mesh->getAABB().getCenter());
        //pos.y= mesh->getAABB().min.y+0.2f;
        pos.x = ((rand()%1024)/1024.0-0.5)*0.25;
        pos.y = mesh->getAABB().min.y+0.2f;
        //pos.y+=((rand()%1024)/1024.0-0.5)*0.5;
        //pos.y+=((rand()%1024)/1024.0-0.5)*0.75;
        //pos.x+=((rand()%1024)/1024.0-0.5)*0.5;
        pos.z=-glm::dvec3(mesh->getAABB().getCenter()).z+((rand()%1024)/1024.0-0.5)*0.25;
        //pos = glm::dvec3(0.0);
        addParticle(Particle(lifeTime*((rand()%1024)/1024.0),pos));
    }
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
    double scale = /*0.5*0.5;*/((decMesh.resolution/2)*(decMesh.resolution/2));
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

            /*
            double s1 = decMesh.getFaceSignum(fit->f1);
            double s2 = decMesh.getFaceSignum(fit->f2);
            double s3 = decMesh.getFaceSignum(fit->f3);
            double s4 = decMesh.getFaceSignum(fit->f4);
            double s5 = decMesh.getFaceSignum(fit->f5);
            double s6 = decMesh.getFaceSignum(fit->f6);*/


            double s1 = 1;
            double s2 = 1;
            double s3 = 1;
            double s4 = 1;
            double s5 = 1;
            double s6 = 1;

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

                    /*
                    wedges[i](ie1,j) += (0.25)*(vel1a*vel4b-vel1b*vel4a)*(glm::dot(glm::cross(-f1.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie2,j) += (0.25)*(vel1a*vel6b-vel1b*vel6a)*(glm::dot(glm::cross(f1.normal,-f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (0.25)*(vel1a*vel3b-vel1b*vel3a)*(glm::dot(glm::cross(f1.normal,-f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (0.25)*(vel1a*vel5b-vel1b*vel5a)*(glm::dot(glm::cross(-f1.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (0.25)*(vel2a*vel4b-vel2b*vel4a)*(glm::dot(glm::cross(-f2.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (0.25)*(vel2a*vel5b-vel2b*vel5a)*(glm::dot(glm::cross(f2.normal,-f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (0.25)*(vel2a*vel3b-vel2b*vel3a)*(glm::dot(glm::cross(f2.normal,-f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (0.25)*(vel2a*vel6b-vel2b*vel6a)*(glm::dot(glm::cross(-f2.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (0.25)*(vel5a*vel4b-vel5b*vel4a)*(glm::dot(glm::cross(-f5.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (0.25)*(vel5a*vel3b-vel5b*vel3a)*(glm::dot(glm::cross(f5.normal,-f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (0.25)*(vel6a*vel4b-vel6b*vel4a)*(glm::dot(glm::cross(-f6.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (0.25)*(vel6a*vel3b-vel6b*vel3a)*(glm::dot(glm::cross(f6.normal,-f3.normal),glm::dvec3(1.0,1.0,1.0)));*/


                    wedges[i](ie1,j) += (scale)*(vel1a*vel4b-vel1b*vel4a)*(glm::dot(glm::cross(f1.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie2,j) += (scale)*(vel6a*vel1b-vel6b*vel1a)*(glm::dot(glm::cross(f6.normal,f1.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (scale)*(vel3a*vel1b-vel3b*vel1a)*(glm::dot(glm::cross(f3.normal,f1.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (scale)*(vel1a*vel5b-vel1b*vel5a)*(glm::dot(glm::cross(f1.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (scale)*(vel4a*vel2b-vel4b*vel2a)*(glm::dot(glm::cross(f4.normal,f2.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (scale)*(vel5a*vel2b-vel5b*vel2a)*(glm::dot(glm::cross(f5.normal,f2.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (scale)*(vel2a*vel3b-vel2b*vel3a)*(glm::dot(glm::cross(f2.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (scale)*(vel2a*vel6b-vel2b*vel6a)*(glm::dot(glm::cross(f2.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (scale)*(vel5a*vel4b-vel5b*vel4a)*(glm::dot(glm::cross(f5.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (scale)*(vel3a*vel5b-vel3b*vel5a)*(glm::dot(glm::cross(f3.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (scale)*(vel4a*vel6b-vel4b*vel6a)*(glm::dot(glm::cross(f4.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (scale)*(vel6a*vel3b-vel6b*vel3a)*(glm::dot(glm::cross(f6.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                }
            }
        }
    }

    for(unsigned int i=0;i<nEigenFunctions;i++)
    {
        for(unsigned int j=0;j<nEigenFunctions;j++)
        {
            advection[i].col(j) = (1.0/eigenValues(i))*wedges[i].col(j);
        }
    }

    for(unsigned int i=0;i<eigenValues.rows();i++)
    {
        advection[i] = vortBasisField.transpose()*advection[i];
    }
}

void PSFSolver::drawParticles(ShaderProgram* program,const glm::mat4& pvm)
{
    particles->syncGPU();
    Particle::setVertexAttribs();
    Particle::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    program->uploadVec4("color",glm::vec4(0.0,0.1,0.0,1.0));
    program->uploadScalar("lifeTime",lifeTime);
    glDrawArrays(GL_POINTS,0,particles->getNumParticles());
}
