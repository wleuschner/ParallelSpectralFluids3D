#include"psfsolvergpu.h"
#include <GL/glew.h>
#include <viennacl/linalg/lanczos.hpp>
#include <viennacl/scalar.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/ocl/backend.hpp>
#include <fstream>
#include <string>
#include <chrono>
#include <CL/cl_gl.h>
#include <CL/cl_gl_ext.h>
#include <GL/glx.h>
#include <QElapsedTimer>
#include <Eigen/Core>
#undef Complex
#undef Success
#include "../../Spectra/MatOp/SparseGenRealShiftSolve.h"
#include "../../Spectra/GenEigsRealShiftSolver.h"
#include "../../DEC/dec.h"
#define WARP_SHIFT 4
#define GRP_SHIFT 4
#define BANK_OFFSET(n)     (((n) >> WARP_SHIFT) + ((n) >> GRP_SHIFT))

PSFSolverGPU::PSFSolverGPU(cl_context context,cl_device_id device,cl_command_queue queue) : AbstractSolver()
{
    cl_context_id = context;
    cl_queue = queue;
    device_id = device;
    viennacl::ocl::setup_context(0,cl_context_id,device_id,cl_queue);

    std::ifstream f("Res/Simulation/simulation.cl");
    std::string source((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    cl_int ret;

    const char* src = source.data();
    program = clCreateProgramWithSource(cl_context_id,1,&src,0,&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not Create Program"<<std::endl;
        exit(-1);
    }

    ret = clBuildProgram(program,1,&device_id,"-cl-fp32-correctly-rounded-divide-sqrt -cl-no-signed-zeros",0,0);

    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build program"<<std::endl;
        char buffer[2048];
        size_t length;
        clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,sizeof(buffer),buffer,&length);
        std::cout<<buffer<<std::endl;
        exit(-1);
    }

    interp_kernel = clCreateKernel(program,"advection",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    visc_kernel = clCreateKernel(program,"normalization_viscocity_gravity",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    vel_update_kernel = clCreateKernel(program,"update_vel",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    advection_reduce_x_kernel = clCreateKernel(program,"advection_reduce_x",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    advection_reduce_y_kernel = clCreateKernel(program,"advection_reduce_y",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    vel_field_reconstruct_kernel = clCreateKernel(program,"reconstruct_velocity_field",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    energy_kernel = clCreateKernel(program,"energy",&ret);
    if(ret!=CL_SUCCESS)
    {
        std::cout<<"Could not build kernel"<<std::endl;
        exit(-1);
    }

    //Load History Compute Shader
    Shader hist_compute(GL_COMPUTE_SHADER,"Res/Volume/history.comp");
    if(!hist_compute.compile())
    {
        std::cout<<hist_compute.compileLog()<<std::endl;
        exit(-1);
    }
    std::cout<<hist_compute.compileLog()<<std::endl;
    historyComputeShader = new ShaderProgram();
    historyComputeShader->attachShader(hist_compute);
    if(!historyComputeShader->link())
    {
        std::cout<<historyComputeShader->linkLog()<<std::endl;
        exit(-1);
    }

    //Load Volume Compute Shader
    Shader vol_compute(GL_COMPUTE_SHADER,"Res/Volume/volume.comp");
    if(!vol_compute.compile())
    {
        std::cout<<vol_compute.compileLog()<<std::endl;
        exit(-1);
    }
    std::cout<<vol_compute.compileLog()<<std::endl;

    volumeComputeShader = new ShaderProgram();
    volumeComputeShader->attachShader(vol_compute);
    if(!volumeComputeShader->link())
    {
        std::cout<<volumeComputeShader->linkLog()<<std::endl;
        exit(-1);
    }
    std::cout<<volumeComputeShader->linkLog()<<std::endl;

    //Load Volume Blur Compute Shader
    Shader vol_blurx_compute(GL_COMPUTE_SHADER,"Res/Volume/blur_x.comp");
    if(!vol_blurx_compute.compile())
    {
        std::cout<<vol_blurx_compute.compileLog()<<std::endl;
        exit(-1);
    }
    std::cout<<vol_blurx_compute.compileLog()<<std::endl;

    volumeBlurXComputeShader = new ShaderProgram();
    volumeBlurXComputeShader->attachShader(vol_compute);
    if(!volumeBlurXComputeShader->link())
    {
        std::cout<<volumeBlurXComputeShader->linkLog()<<std::endl;
        exit(-1);
    }
    std::cout<<volumeBlurXComputeShader->linkLog()<<std::endl;

    Shader vol_blury_compute(GL_COMPUTE_SHADER,"Res/Volume/blur_y.comp");
    if(!vol_blury_compute.compile())
    {
        std::cout<<vol_blury_compute.compileLog()<<std::endl;
        exit(-1);
    }
    std::cout<<vol_blury_compute.compileLog()<<std::endl;

    volumeBlurYComputeShader = new ShaderProgram();
    volumeBlurYComputeShader->attachShader(vol_compute);
    if(!volumeBlurYComputeShader->link())
    {
        std::cout<<volumeBlurYComputeShader->linkLog()<<std::endl;
        exit(-1);
    }
    std::cout<<volumeBlurYComputeShader->linkLog()<<std::endl;

    Shader vol_blurz_compute(GL_COMPUTE_SHADER,"Res/Volume/blur_z.comp");
    if(!vol_blurx_compute.compile())
    {
        std::cout<<vol_blurz_compute.compileLog()<<std::endl;
        exit(-1);
    }
    std::cout<<vol_blurz_compute.compileLog()<<std::endl;

    volumeBlurZComputeShader = new ShaderProgram();
    volumeBlurZComputeShader->attachShader(vol_compute);
    if(!volumeBlurZComputeShader->link())
    {
        std::cout<<volumeBlurZComputeShader->linkLog()<<std::endl;
        exit(-1);
    }
    std::cout<<volumeBlurZComputeShader->linkLog()<<std::endl;


    //Load Blit Shader
    Shader blitVert(GL_VERTEX_SHADER,"Res/Effects/Volume/volume.vert");
    if(!blitVert.compile())
    {
        std::cout<<blitVert.compileLog()<<std::endl;
    }
    Shader blitFrag(GL_FRAGMENT_SHADER,"Res/Effects/Volume/blit.frag");
    if(!blitFrag.compile())
    {
        std::cout<<blitFrag.compileLog()<<std::endl;
    }
    blitProgram = new ShaderProgram();
    blitProgram->attachShader(blitVert);
    blitProgram->attachShader(blitFrag);
    if(!blitProgram->link())
    {
        std::cout<<blitProgram->linkLog()<<std::endl;
    }
    blitProgram->bind();

    refreshParticles = 60*10;
}

void PSFSolverGPU::integrate()
{
    glFinish();
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
        calculateEnergy(0);
        calculateVelocity();
        calculateEnergy(1);
        updateVelocity();
        externalForces();
        reconstructVelocityField();
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
        updateParticles();
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        if(isBenchmark)
        {
            double t = std::chrono::duration<double,std::milli>(end-begin).count();
            benchmarkSums[1] +=t;
            benchmark_file<<t<<std::endl;
        }
    }
    clFinish(cl_queue);
}

void PSFSolverGPU::calculateEnergy(unsigned int index)
{
    cl_int res;
    size_t global_work_size[] = {nEigenFunctionsAligned/8,1,1};
    size_t local_work_size[] = {nEigenFunctionsAligned/8,1,1};

    clSetKernelArg(energy_kernel,0,sizeof(cl_mem),&basisCoeffBuffer);
    clSetKernelArg(energy_kernel,1,sizeof(cl_mem),&energyBuffer);
    clSetKernelArg(energy_kernel,2,sizeof(cl_uint),&index);
    clSetKernelArg(energy_kernel,3,sizeof(cl_double)*(2*local_work_size[0]+BANK_OFFSET(2*local_work_size[0])),0);
    res = clEnqueueNDRangeKernel(cl_queue,energy_kernel,1,0,global_work_size,local_work_size,0,0,0);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }
}

void PSFSolverGPU::calculateVelocity()
{
    cl_event event;
    cl_int res;

    size_t advection_reduce_work_size[] = {nEigenFunctionsAligned/8u,
                                           nEigenFunctions,
                                           nEigenFunctions};
    size_t advection_reduce_local_size[] = {nEigenFunctionsAligned/8u,1,1};

    clSetKernelArg(advection_reduce_x_kernel,0,sizeof(cl_mem),&advectionMatrices);
    clSetKernelArg(advection_reduce_x_kernel,1,sizeof(cl_mem),&advectionScratchBuffer);
    clSetKernelArg(advection_reduce_x_kernel,2,sizeof(cl_mem),&basisCoeffBuffer);
    clSetKernelArg(advection_reduce_x_kernel,3,sizeof(cl_double)*(2*advection_reduce_local_size[0]+BANK_OFFSET(2*advection_reduce_local_size[0])),NULL);
    res = clEnqueueNDRangeKernel(cl_queue,advection_reduce_x_kernel,3,0,advection_reduce_work_size,advection_reduce_local_size,0,0,&event);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }

    clSetKernelArg(advection_reduce_y_kernel,0,sizeof(cl_mem),&advectionScratchBuffer);
    clSetKernelArg(advection_reduce_y_kernel,1,sizeof(cl_mem),&velocityUpdateBuffer);
    clSetKernelArg(advection_reduce_y_kernel,2,sizeof(cl_mem),&basisCoeffBuffer);
    clSetKernelArg(advection_reduce_y_kernel,3,sizeof(cl_double)*(2*advection_reduce_local_size[0]+BANK_OFFSET(2*advection_reduce_local_size[0])),NULL);
    res = clEnqueueNDRangeKernel(cl_queue,advection_reduce_y_kernel,2,0,advection_reduce_work_size,advection_reduce_local_size,0,0,0);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }
}

void PSFSolverGPU::updateVelocity()
{
    cl_int res;
    size_t global_work_size = nEigenFunctionsAligned/4;
    clSetKernelArg(vel_update_kernel,0,sizeof(cl_mem),&basisCoeffBuffer);
    clSetKernelArg(vel_update_kernel,1,sizeof(cl_mem),&velocityUpdateBuffer);
    clSetKernelArg(vel_update_kernel,2,sizeof(cl_double),&timeStep);
    res = clEnqueueNDRangeKernel(cl_queue,vel_update_kernel,1,0,&global_work_size,0,0,0,0);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }
}


void PSFSolverGPU::externalForces()
{
    cl_int res;
    size_t global_work_size = nEigenFunctionsAligned/4;
    double argGravity = gravityActive==true?1.0:0.0;

    clSetKernelArg(visc_kernel,0,sizeof(cl_mem),&energyBuffer);
    clSetKernelArg(visc_kernel,1,sizeof(cl_mem),&basisCoeffBuffer);
    clSetKernelArg(visc_kernel,2,sizeof(cl_mem),&eigenValuesBuffer);
    clSetKernelArg(visc_kernel,3,sizeof(cl_mem),&gravityBuffer);
    clSetKernelArg(visc_kernel,4,sizeof(cl_double),&viscosity);
    clSetKernelArg(visc_kernel,5,sizeof(cl_double),&timeStep);
    clSetKernelArg(visc_kernel,6,sizeof(cl_double),&argGravity);
    res = clEnqueueNDRangeKernel(cl_queue,visc_kernel,1,0,&global_work_size,0,0,0,0);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }
}

void PSFSolverGPU::reconstructVelocityField()
{
    cl_int res;
    size_t faces = velBasisField.rows();
    size_t global_work_size[] = {nEigenFunctionsAligned/8,faces};
    size_t local_work_size[] = {nEigenFunctionsAligned/8,1,1};
    clSetKernelArg(vel_field_reconstruct_kernel,0,sizeof(cl_mem),&velocityBasisFieldBuffer);
    clSetKernelArg(vel_field_reconstruct_kernel,1,sizeof(cl_mem),&velocityFieldBuffer);
    clSetKernelArg(vel_field_reconstruct_kernel,2,sizeof(cl_mem),&basisCoeffBuffer);
    clSetKernelArg(vel_field_reconstruct_kernel,3,sizeof(cl_double)*(2*local_work_size[0]+BANK_OFFSET(2*local_work_size[0])),NULL);
    res = clEnqueueNDRangeKernel(cl_queue,vel_field_reconstruct_kernel,2,0,global_work_size,local_work_size,0,0,0);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }
}

void PSFSolverGPU::updateParticles()
{
    cl_int res;
    cl_event event;

    particlesBuffer = clCreateFromGLBuffer(cl_context_id,CL_MEM_READ_WRITE,particles->id,&res);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Could not create CLGL Buffer Object!"<<std::endl;
        exit(-1);
    }

    res = clEnqueueAcquireGLObjects(cl_queue,1,&particlesBuffer,0,0,&event);
    res = clWaitForEvents(1,&event);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Failed to Acquire GL objects";
        exit(-1);
    }
    glm::uvec3 random = glm::uvec3(rand(),rand(),rand());
    glm::uvec4 dims = glm::uvec4(decMesh.getDimensions(),0);
    glm::vec4 aabb_min = glm::vec4(mesh->getAABB().min,0.0);
    glm::vec4 aabb_max = glm::vec4(mesh->getAABB().max,0.0);
    glm::vec4 aabb_extent = glm::vec4(mesh->getAABB().getExtent(),0.0);
    float resf = resolution;
    cl_float timestep = timeStep;
    size_t global_work_size = maxParticles;
    clSetKernelArg(interp_kernel,0,sizeof(cl_mem),&particlesBuffer);
    clSetKernelArg(interp_kernel,1,sizeof(cl_mem),&velocityFieldBuffer);
    clSetKernelArg(interp_kernel,2,sizeof(cl_mem),&signBitStringHandle);
    clSetKernelArg(interp_kernel,3,sizeof(cl_uint4),&dims);
    clSetKernelArg(interp_kernel,4,sizeof(cl_float4),&aabb_min);
    clSetKernelArg(interp_kernel,5,sizeof(cl_float4),&aabb_max);
    clSetKernelArg(interp_kernel,6,sizeof(cl_float4),&aabb_extent);
    clSetKernelArg(interp_kernel,7,sizeof(cl_float),&resf);
    clSetKernelArg(interp_kernel,8,sizeof(cl_float),&timestep);
    clSetKernelArg(interp_kernel,9,sizeof(cl_uint3),&random);
    clSetKernelArg(interp_kernel,10,sizeof(cl_float),&lifeTime);
    res = clEnqueueNDRangeKernel(cl_queue,interp_kernel,1,0,&global_work_size,0,0,0,&event);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel"<<std::endl;
    }
    res = clWaitForEvents(1,&event);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Could not wait for kernel completion"<<std::endl;
    }
    res = clEnqueueReleaseGLObjects(cl_queue,1,&particlesBuffer,0,0,&event);
    res = clWaitForEvents(1,&event);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Failed to Release GL objects";
        exit(-1);
    }
    clReleaseMemObject(particlesBuffer);
}

void PSFSolverGPU::buildLaplace()
{
    Eigen::SparseMatrix<double> mat = 1.0*derivative1(decMesh,false)*hodge2(decMesh,true)*derivative1(decMesh,true)*hodge2(decMesh,false);
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
                    if(backProjection.isApprox(tempEigenValues(i).real()*tempEigenVectors.col(i).real(),0.1))
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
                    gravity(it->id) = (it->normal.y>0?1:-1)*9.81;
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

                    vorticityField(it->id) = ((p1-p2).z>0?1.0:-1.0)*100;
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
                   it->center.y<aabb.min.y+0.8 /*&& abs(it->center.x)<0.4 && abs(it->center.z)<0.4*/)
                        velocityField(decMesh.getFaceIndex(*it)) = (it->normal.y>0?-1:1)*1.0;
            }
        }
    }
    setInitialVelocityField(velocityField);

    /*for(unsigned int i=0;i<400000;i++)
    {
        glm::dvec3 pos = glm::dvec3(mesh->getAABB().getCenter());
        //pos.y= mesh->getAABB().min.y+0.2f;
        pos.y+=((rand()%1024)/1024.0-0.5)*0.25;
        //pos.y+=((rand()%1024)/1024.0-0.5)*0.75;
        //pos.x+=((rand()%1024)/1024.0-0.5)*1.0;
        pos.x+=((rand()%1024)/1024.0-0.5)*0.25;
        pos.z=-glm::dvec3(mesh->getAABB().getCenter()).z+((rand()%1024)/1024.0-0.5)*0.25;
        //pos.z=-glm::dvec3(mesh->getAABB().getCenter()).z+((rand()%1024)/1024.0-0.5)*1.0;
        //pos = glm::dvec3(0.0);
        addParticle(Particle(lifeTime*((rand()%1024)/1024.0),pos));
    }*/
}

void PSFSolverGPU::buildAdvection()
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
    uint nextPower2 = pow(2, ceil(log(nEigenFunctions)/log(2)));
    nEigenFunctionsAligned = nEigenFunctions + (nextPower2-nEigenFunctions%nextPower2)*(nEigenFunctions%nextPower2==0?0:1);
    uploadFaceSigns();
    uploadAdvectionMatrices();
    uploadVelocityBasisField();
    uploadVelocityField();
    uploadBasisCoeffs();
    uploadGravity();
    uploadEigenValues();
    uploadEnergy();
    particles->syncGPU();
}

void PSFSolverGPU::uploadAdvectionMatrices()
{
    unsigned int gpuUploadSize = nEigenFunctions*nEigenFunctions*nEigenFunctionsAligned;
    std::vector<double> gpuUpload;
    gpuUpload.resize(gpuUploadSize);
    memset(gpuUpload.data(),0,sizeof(double)*gpuUploadSize);
    cl_int ret;
    velocityUpdateBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,(nEigenFunctionsAligned)*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }

    advectionScratchBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,(nEigenFunctions*nEigenFunctionsAligned)*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }


    for(unsigned int i=0;i<nEigenFunctions;i++)
    {
        Eigen::MatrixXd advectionTemp = advection[i].transpose();
        #pragma omp parallel for
        for(unsigned int j=0;j<nEigenFunctions;j++)
        {
            memcpy(gpuUpload.data() + i*(nEigenFunctions*nEigenFunctionsAligned)+j*(nEigenFunctionsAligned),advectionTemp.col(j).data(),sizeof(double)*nEigenFunctions);
        }
    }

    advectionMatrices = clCreateBuffer(cl_context_id,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,gpuUpload.size()*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadBasisCoeffs()
{
    unsigned int gpuUploadSize = nEigenFunctionsAligned;
    std::vector<double> gpuUpload;
    gpuUpload.resize(gpuUploadSize);
    memset(gpuUpload.data(),0,sizeof(double)*gpuUploadSize);
    cl_int ret;
    memcpy(gpuUpload.data(),basisCoeff.data(),sizeof(double)*nEigenFunctions);
    basisCoeffBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,gpuUpload.size()*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadFaceSigns()
{
    std::vector<unsigned int>& signBitString = decMesh.getSignBitString();
    cl_int ret;
    signBitStringHandle = clCreateBuffer(cl_context_id,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,signBitString.size()*sizeof(unsigned int),signBitString.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadVelocityBasisField()
{
    unsigned int gpuUploadSize = velBasisField.rows()*nEigenFunctionsAligned;
    std::vector<double> gpuUpload;
    gpuUpload.resize(gpuUploadSize);
    memset(gpuUpload.data(),0,sizeof(double)*gpuUploadSize);
    cl_int ret;

    Eigen::MatrixXd velBasisFieldTemp = velBasisField.transpose();
    for(unsigned int i=0;i<velBasisFieldTemp.cols();i++)
    {
        memcpy(gpuUpload.data() + i*(nEigenFunctionsAligned),velBasisFieldTemp.col(i).data(),sizeof(double)*nEigenFunctions);
    }

    velocityBasisFieldBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,gpuUpload.size()*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadGravity()
{
    unsigned int gpuUploadSize = nEigenFunctionsAligned;
    std::vector<double> gpuUpload;
    gpuUpload.resize(gpuUploadSize);
    memset(gpuUpload.data(),0,sizeof(double)*gpuUploadSize);
    cl_int ret;
    memcpy(gpuUpload.data(),gravity.data(),sizeof(double)*nEigenFunctions);
    gravityBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,gpuUpload.size()*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadEigenValues()
{
    unsigned int gpuUploadSize = nEigenFunctionsAligned;
    std::vector<double> gpuUpload;
    gpuUpload.resize(gpuUploadSize);
    memset(gpuUpload.data(),0,sizeof(double)*gpuUploadSize);
    cl_int ret;
    memcpy(gpuUpload.data(),eigenValues.data(),sizeof(double)*nEigenFunctions);
    eigenValuesBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,gpuUpload.size()*sizeof(double),gpuUpload.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadVelocityField()
{
    unsigned int gpuUploadSize = velBasisField.rows();
    cl_int ret;
    velocityFieldBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_WRITE,gpuUploadSize*sizeof(float),0,&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::uploadEnergy()
{
    cl_int ret;
    energyBuffer = clCreateBuffer(cl_context_id,CL_MEM_READ_WRITE,2*sizeof(double),0,&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }
}

void PSFSolverGPU::drawParticles(ShaderProgram* program,const glm::mat4& pvm)
{
    /*
    particles->bind();
    Particle::setVertexAttribs();
    Particle::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    program->uploadVec4("color",glm::vec4(0.0,0.1,0.0,1.0));
    glDrawArrays(GL_POINTS,0,particles->getNumParticles());
    return;*/

    float rnd = (rand()%(16<<2))/float(16<<2);
    glm::vec3 jitter = glm::normalize(glm::vec3((rand()%1024)/1024.0,(rand()%1024)/1024.0,(rand()%1024)/1024.0));
    glm::vec3 dims2 = glm::vec3(128,128,128);
    glm::vec3 extent = mesh->getAABB().getExtent();

    histogramTexture->clearImage();
    particles->bindCompute(0);
    histogramTexture->bindCompute(1);
    historyComputeShader->bind();
    historyComputeShader->uploadVec3("aabb_min",mesh->getAABB().min);
    historyComputeShader->uploadVec3("aabb_extent",extent);
    historyComputeShader->uploadInt("volumeTexture",1);
    historyComputeShader->dispatch(maxParticles/1024,1,1,maxParticles,1,1);
    GLsync syncObj;
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    syncObj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE,0);
    glClientWaitSync(syncObj,0,1000*1000*1000*2);
    glDeleteSync(syncObj);

    histogramTexture->bindCompute(0);
    volumeTextures[0]->bindFloatCompute(1);
    volumeComputeShader->bind();
    volumeComputeShader->uploadInt("histogramTexture",0);
    volumeComputeShader->uploadInt("volumeTexture",1);
    volumeComputeShader->uploadScalar("lifeTime",lifeTime);
    volumeComputeShader->dispatch(256,256,1,1,1,1);

    int result;
    if(result!=GL_NO_ERROR)
    {
        std::cout<<"Compute Error"<<result<<std::endl;
    }
    glMemoryBarrier(GL_ALL_BARRIER_BITS);
    syncObj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE,0);
    glClientWaitSync(syncObj,0,1000*1000*1000*2);
    glDeleteSync(syncObj);

    unsigned int currentOutTex = 0;
    /*for(unsigned int i=0;i<8;i++)
    {
        volumeTextures[currentOutTex]->bindFloatCompute(0);
        volumeTextures[currentOutTex^1]->bindFloatCompute(1);
        volumeBlurXComputeShader->bind();
        volumeBlurXComputeShader->uploadInt("histogramTexture",0);
        volumeBlurXComputeShader->uploadInt("volumeTexture",1);
        volumeBlurXComputeShader->dispatch(1024,1024,1,1,1,1);

        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        syncObj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE,0);
        glClientWaitSync(syncObj,0,1000*1000*1000*2);
        glDeleteSync(syncObj);
        currentOutTex^=1;

        volumeTextures[currentOutTex]->bindFloatCompute(0);
        volumeTextures[currentOutTex^1]->bindFloatCompute(1);
        volumeBlurYComputeShader->bind();
        volumeBlurYComputeShader->uploadInt("histogramTexture",0);
        volumeBlurYComputeShader->uploadInt("volumeTexture",1);
        volumeBlurYComputeShader->dispatch(1024,1024,1,1,1,1);

        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        syncObj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE,0);
        glClientWaitSync(syncObj,0,1000*1000*1000*2);
        glDeleteSync(syncObj);
        currentOutTex^=1;

        volumeTextures[currentOutTex]->bindFloatCompute(0);
        volumeTextures[currentOutTex^1]->bindFloatCompute(1);
        volumeBlurZComputeShader->bind();
        volumeBlurZComputeShader->uploadInt("histogramTexture",0);
        volumeBlurZComputeShader->uploadInt("volumeTexture",1);
        volumeBlurZComputeShader->dispatch(1024,1024,1,1,1,1);

        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        syncObj = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE,0);
        glClientWaitSync(syncObj,0,1000*1000*1000*2);
        glDeleteSync(syncObj);
        currentOutTex^=1;
    }*/

    historyBuffer[currentHistory]->bind();
    glClear(GL_COLOR_BUFFER_BIT);
    glBindBuffer(GL_ARRAY_BUFFER,fullscreenVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fullscreenIBO);
    glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,sizeof(glm::vec4),(void*)0);
    volumeTextures[0]->bind(0);
    program->bind();
    program->uploadScalar("rnd",rnd);
    program->uploadVec3("jitter",jitter);
    program->uploadVec4("viewport_size",viewport_size);
    program->uploadVec3("camera_position",camera_position);
    program->uploadMat4("view_mat",view_mat);
    program->uploadScalar("step_size",extent.x/256.0);
    program->uploadVec3("aabb_min",mesh->getAABB().min);
    program->uploadVec3("aabb_max",mesh->getAABB().max);
    program->uploadInt("volumeTexture",0);
    program->uploadMat4("pvm",pvm);
    program->uploadLight("light",light,view_mat);
    glDrawElements(GL_TRIANGLES,6,GL_UNSIGNED_INT,0);


    historyBuffer[currentHistory]->unbind();
    glBindBuffer(GL_ARRAY_BUFFER,fullscreenVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fullscreenIBO);
    glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE,sizeof(glm::vec4),(void*)0);
    colorAttachments[0]->bind(0);
    colorAttachments[1]->bind(1);
    colorAttachments[2]->bind(2);
    colorAttachments[3]->bind(3);
    colorAttachments[4]->bind(4);
    colorAttachments[5]->bind(5);
    colorAttachments[6]->bind(6);
    colorAttachments[7]->bind(7);
    blitProgram->bind();
    blitProgram->uploadVec4("viewport_size",viewport_size);
    blitProgram->uploadInt("texture_a",(currentHistory)%8);
    blitProgram->uploadInt("texture_b",(currentHistory+1)%8);
    blitProgram->uploadInt("texture_c",(currentHistory+2)%8);
    blitProgram->uploadInt("texture_d",(currentHistory+3)%8);
    blitProgram->uploadInt("texture_e",(currentHistory+4)%8);
    blitProgram->uploadInt("texture_f",(currentHistory+5)%8);
    blitProgram->uploadInt("texture_g",(currentHistory+6)%8);
    blitProgram->uploadInt("texture_h",(currentHistory+7)%8);
    glDrawElements(GL_TRIANGLES,6,GL_UNSIGNED_INT,0);
    currentHistory=(currentHistory+1)%8;

    /*particles->bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    program->uploadVec4("color",glm::vec4(0.0,0.1,0.0,1.0));
    program->uploadScalar("lifeTime",lifeTime);
    glDrawArrays(GL_POINTS,0,maxParticles);*/
}
