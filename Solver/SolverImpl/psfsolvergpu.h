#ifndef __PSF_SOLVER_GPU_H_
#define __PSF_SOLVER_GPU_H_
#include "../abstractsolver.h"
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/ocl/backend.hpp>
#include <CL/cl.h>

class PSFSolverGPU : public AbstractSolver
{
public:
    PSFSolverGPU(cl_context context,cl_device_id device,cl_command_queue queue);
    void integrate();

    void drawParticles(ShaderProgram* program,const glm::mat4& pvm);
    Light light;

protected:
    void buildLaplace();
    void buildAdvection();
private:
    void calculateEnergy(unsigned int index);
    void calculateVelocity();
    void updateVelocity();
    void externalForces();
    void reconstructVelocityField();
    void updateParticles();


    void uploadAdvectionMatrices();
    void uploadBasisCoeffs();
    void uploadFaceSigns();
    void uploadVelocityBasisField();
    void uploadGravity();
    void uploadEigenValues();
    void uploadVelocityField();
    void uploadEnergy();

    ShaderProgram* historyComputeShader;
    ShaderProgram* volumeComputeShader;
    ShaderProgram* volumeBlurXComputeShader;
    ShaderProgram* volumeBlurYComputeShader;
    ShaderProgram* volumeBlurZComputeShader;
    ShaderProgram* blitProgram;

    unsigned int nEigenFunctionsAligned;

    uint32_t refreshParticles;
    cl_context cl_context_id;
    cl_command_queue cl_queue;
    cl_device_id device_id;
    cl_program program;
    cl_kernel interp_kernel;
    cl_kernel visc_kernel;
    cl_kernel vel_update_kernel;
    cl_kernel advection_reduce_x_kernel;
    cl_kernel advection_reduce_y_kernel;
    cl_kernel vel_field_reconstruct_kernel;
    cl_kernel energy_kernel;

    cl_mem signBitStringHandle;
    cl_mem advectionMatrices;
    cl_mem advectionScratchBuffer;
    cl_mem velocityUpdateBuffer;
    cl_mem velocityFieldBuffer;
    cl_mem basisCoeffBuffer;
    cl_mem velocityBasisFieldBuffer;
    cl_mem eigenValuesBuffer;
    cl_mem gravityBuffer;
    cl_mem energyBuffer;

    cl_mem particlesBuffer;

    viennacl::compressed_matrix<double> vcl_curl;

    viennacl::vector<double> vcl_vorticityField;
    viennacl::vector<double> vcl_velocityField;

    viennacl::scalar<double> vcl_e1;
    viennacl::scalar<double> vcl_e2;
    viennacl::vector<double> vcl_velocity;
    viennacl::vector<double> vcl_gravity;
    viennacl::vector<double> vcl_eigenValues;
    viennacl::vector<double> vcl_basisCoeff;
    viennacl::matrix<double> vcl_velBasisField;
    viennacl::matrix<double> vcl_vortBasisField;
};

#endif
