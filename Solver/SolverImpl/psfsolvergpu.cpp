#include"psfsolvergpu.h"
#include <GL/glew.h>
#include <viennacl/linalg/lanczos.hpp>
#include <viennacl/scalar.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/ocl/backend.hpp>
#include <fstream>
#include <string>
#include <CL/cl_gl.h>
#include <CL/cl_gl_ext.h>
#include <GL/glx.h>
#include <QElapsedTimer>
#include "../../Spectra/MatOp/SparseSymShiftSolve.h"
#include "../../Spectra/SymEigsShiftSolver.h"
#include "../../DEC/dec.h"

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
    refreshParticles = 60*10;
}

void PSFSolverGPU::integrate()
{
    QElapsedTimer timer;
    timer.start();

    cl_int res;
    size_t global_work_size = 1;
    size_t local_work_size = 256;
    cl_mem basisCoeff_handle = vcl_basisCoeff.handle().opencl_handle();
    cl_mem eigenVals_handle = vcl_eigenValues.handle().opencl_handle();
    cl_mem gravity_handle = vcl_gravity.handle().opencl_handle();
    cl_mem vel_handle = vcl_velocityField.handle().opencl_handle();
    cl_mem e1_handle = vcl_e1.handle().opencl_handle();
    cl_mem e2_handle = vcl_e2.handle().opencl_handle();
    double argGravity = gravityActive==true?1.0:0.0;

    vcl_e1 = viennacl::linalg::inner_prod(vcl_basisCoeff,vcl_basisCoeff);

    viennacl::scalar<double> val1 = 0.0;
    cl_mem val1_handle = val1.handle().opencl_handle();
    clSetKernelArg(vel_update_kernel,0,sizeof(cl_mem),&basisCoeff_handle);
    clSetKernelArg(vel_update_kernel,1,sizeof(cl_mem),&val1_handle);
    clSetKernelArg(vel_update_kernel,2,sizeof(cl_double),&timeStep);
    for(unsigned int k=0;k<nEigenFunctions;k++)
    {
        val1 = viennacl::linalg::inner_prod(vcl_basisCoeff,viennacl::linalg::prod(vcl_advection[k],vcl_basisCoeff));
        clSetKernelArg(vel_update_kernel,3,sizeof(cl_uint),&k);
        res = clEnqueueNDRangeKernel(cl_queue,vel_update_kernel,1,0,&global_work_size,0,0,0,0);
        if(res!=CL_SUCCESS)
        {
            std::cout<<"Unable to execute Kernel "<<res<<std::endl;
        }
        //vcl_velocity(k) += viennacl::linalg::inner_prod(vcl_basisCoeff,viennacl::linalg::prod(vcl_advection[k],vcl_basisCoeff));
    }
    //vcl_basisCoeff += vcl_timestep*vcl_velocity;

    vcl_e2 = viennacl::linalg::inner_prod(vcl_basisCoeff,vcl_basisCoeff);

    //basisCoeff *= vcl_e1/vcl_e2;

    global_work_size = nEigenFunctions;
    clSetKernelArg(visc_kernel,0,sizeof(cl_mem),&e1_handle);
    clSetKernelArg(visc_kernel,1,sizeof(cl_mem),&e2_handle);
    clSetKernelArg(visc_kernel,2,sizeof(cl_mem),&basisCoeff_handle);
    clSetKernelArg(visc_kernel,3,sizeof(cl_mem),&eigenVals_handle);
    clSetKernelArg(visc_kernel,4,sizeof(cl_mem),&gravity_handle);
    clSetKernelArg(visc_kernel,5,sizeof(cl_double),&viscosity);
    clSetKernelArg(visc_kernel,6,sizeof(cl_double),&timeStep);
    clSetKernelArg(visc_kernel,7,sizeof(cl_double),&argGravity);
    res = clEnqueueNDRangeKernel(cl_queue,visc_kernel,1,0,&global_work_size,0,0,0,0);
    if(res!=CL_SUCCESS)
    {
        std::cout<<"Unable to execute Kernel "<<res<<std::endl;
    }
    vcl_velocityField = viennacl::linalg::prod(vcl_velBasisField,vcl_basisCoeff);

    cl_int ret;
    unsigned int workGroupSize = 128;
    unsigned int workGroups = std::ceil(maxParticles/float(workGroupSize));
    glFinish();
    cl_event event;

    particlesBuffer = clCreateFromGLBuffer(cl_context_id,CL_MEM_READ_WRITE,particles->id,&ret);
    if(ret!=CL_SUCCESS)
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
    global_work_size = maxParticles;
    local_work_size = {256};
    clSetKernelArg(interp_kernel,0,sizeof(cl_mem),&particlesBuffer);
    clSetKernelArg(interp_kernel,1,sizeof(cl_mem),&vel_handle);
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
    clFinish(cl_queue);
    std::cout<<timer.elapsed()<<std::endl;
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

            int nconv = solver.compute(1000,1e-10,Spectra::WHICH_LM);
            if(solver.info()==Spectra::SUCCESSFUL)
            {

                Eigen::VectorXd tempEigenValues = solver.eigenvalues().real();
                Eigen::MatrixXd tempEigenVectors = solver.eigenvectors().real();
                bool onlyZero=true;
                for(unsigned int i=0;i<nconv && foundEigenValues<nEigenFunctions;i++)
                {
                    Eigen::VectorXd backProjection = mat*tempEigenVectors.col(i);
                    if(backProjection.isApprox(tempEigenValues(i)*tempEigenVectors.col(i),1))
                    //if(std::abs(1.0/tempEigenValues(i))>1e-10 && std::abs(1.0/tempEigenValues(i))<1e+10)
                    {
                        bool doubleEv = false;
                        for(unsigned int j=0;j<foundEigenValues;j++)
                        {
                            //if(fabs(eigenValues(j)-tempEigenValues(i))<=std::numeric_limits<float>::epsilon())
                            if(velBasisField.col(j).isApprox(tempEigenVectors.col(i)))
                            {
                                std::cout<<"Already Inside "<<tempEigenValues(i)<<std::endl;
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
                        //std::cout<<"ZERO GARBAGE "<<tempEigenValues(i)<<std::endl;
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
                    gravity(it->id) = (it->normal.y>0?1:-1)*9.81;
                }
            }
        }
    }
    gravity = velBasisField.transpose()*gravity;

    std::vector<unsigned int>& signBitString = decMesh.getSignBitString();
    cl_int ret;
    signBitStringHandle = clCreateBuffer(cl_context_id,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,signBitString.size()*sizeof(unsigned int),signBitString.data(),&ret);
    if(ret!=CL_SUCCESS)
    {
        printf("Could not allocate Buffer");
        exit(ret);
    }

    std::vector<Vertex> vertices = mesh->getVertices();

    vorticityField = Eigen::VectorXd::Zero(decMesh.getNumEdges());
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
                if(glm::dot(it->normal,glm::dvec3(0.0,1.0,0.0)) &&
                   it->center.y<aabb.min.y+0.4 && abs(it->center.x)<0.4 && abs(it->center.z)<0.4)
                        velocityField(decMesh.getFaceIndex(*it)) = (it->normal.y>0?-1:1)*50.0;

            }
        }
    }
    setInitialVelocityField(velocityField);
    vcl_vortBasisField.resize(vortBasisField.rows(),vortBasisField.cols());
    vcl_velBasisField.resize(velBasisField.rows(),velBasisField.cols());
    vcl_gravity.resize(gravity.size());
    vcl_vorticityField.resize(vorticityField.rows(),vorticityField.cols());
    vcl_velocityField.resize(velocityField.rows(),velocityField.cols());
    vcl_basisCoeff.resize(basisCoeff.size());
    vcl_eigenValues.resize(eigenValues.size());

    viennacl::copy(vortBasisField,vcl_vortBasisField);
    viennacl::copy(velBasisField,vcl_velBasisField);
    viennacl::copy(gravity,vcl_gravity);
    viennacl::copy(vorticityField,vcl_vorticityField);
    viennacl::copy(velocityField,vcl_velocityField);
    viennacl::copy(eigenValues,vcl_eigenValues);
    viennacl::copy(basisCoeff,vcl_basisCoeff);
    vcl_velocity = viennacl::zero_vector<double>(nEigenFunctions);
    vcl_e1 = 0.0;
    vcl_e2 = 0.0;
    vcl_timestep = timeStep;

    for(unsigned int i=0;i<400000;i++)
    {
        glm::dvec3 pos = glm::dvec3(mesh->getAABB().getCenter());
        pos.y= mesh->getAABB().min.y+0.2f;
        //pos.y+=((rand()%1024)/1024.0-0.5)*0.5;
        //pos.y+=((rand()%1024)/1024.0-0.5)*0.75;
        pos.x+=((rand()%1024)/1024.0-0.5)*1.0;
        pos.z=-glm::dvec3(mesh->getAABB().getCenter()).z+((rand()%1024)/1024.0-0.5)*1.0;
        //pos = glm::dvec3(0.0);
        addParticle(Particle(lifeTime*((rand()%1024)/1024.0),pos));
    }
    particles->syncGPU();
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
    double scale = ((decMesh.resolution/2)*(decMesh.resolution/2));
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
                    wedges[i](ie2,j) += (scale)*(vel1a*vel6b-vel1b*vel6a)*(glm::dot(glm::cross(f1.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie3,j) += (scale)*(vel1a*vel3b-vel1b*vel3a)*(glm::dot(glm::cross(f1.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie4,j) += (scale)*(vel1a*vel5b-vel1b*vel5a)*(glm::dot(glm::cross(f1.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie5,j) += (scale)*(vel2a*vel4b-vel2b*vel4a)*(glm::dot(glm::cross(f2.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie6,j) += (scale)*(vel2a*vel5b-vel2b*vel5a)*(glm::dot(glm::cross(f2.normal,f5.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie7,j) += (scale)*(vel2a*vel3b-vel2b*vel3a)*(glm::dot(glm::cross(f2.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie8,j) += (scale)*(vel2a*vel6b-vel2b*vel6a)*(glm::dot(glm::cross(f2.normal,f6.normal),glm::dvec3(1.0,1.0,1.0)));

                    wedges[i](ie9,j) += (scale)*(vel5a*vel4b-vel5b*vel4a)*(glm::dot(glm::cross(f5.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie10,j) += (scale)*(vel5a*vel3b-vel5b*vel3a)*(glm::dot(glm::cross(f5.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie11,j) += (scale)*(vel6a*vel4b-vel6b*vel4a)*(glm::dot(glm::cross(f6.normal,f4.normal),glm::dvec3(1.0,1.0,1.0)));
                    wedges[i](ie12,j) += (scale)*(vel6a*vel3b-vel6b*vel3a)*(glm::dot(glm::cross(f6.normal,f3.normal),glm::dvec3(1.0,1.0,1.0)));
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

void PSFSolverGPU::drawParticles(ShaderProgram* program,const glm::mat4& pvm)
{
    glEnableClientState(GL_VERTEX_ARRAY);
    particles->bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    program->uploadVec4("color",glm::vec4(0.0,0.1,0.0,1.0));
    program->uploadScalar("lifeTime",lifeTime);
    glDrawArrays(GL_POINTS,0,maxParticles);
}
