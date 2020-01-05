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
protected:
    void buildLaplace();
    void buildAdvection();
private:
    uint32_t refreshParticles;
    cl_context cl_context_id;
    cl_command_queue cl_queue;
    cl_device_id device_id;
    cl_program program;
    cl_kernel interp_kernel;
    cl_kernel visc_kernel;
    cl_kernel vel_update_kernel;

    cl_mem signBitStringHandle;

    viennacl::ocl::program vcl_prog_psf;
    viennacl::ocl::kernel advection_kernel;
    cl_mem particlesBuffer;

    viennacl::compressed_matrix<double> vcl_curl;

    std::vector<viennacl::vector<double>> vcl_eigenFunctions;

    viennacl::vector<double> vcl_vorticityField;
    viennacl::vector<double> vcl_velocityField;

    std::vector<viennacl::matrix<double>> vcl_advection;

    viennacl::scalar<double> vcl_e1;
    viennacl::scalar<double> vcl_e2;
    viennacl::vector<double> vcl_velocity;
    viennacl::scalar<double> vcl_timestep;
    viennacl::vector<double> vcl_gravity;
    viennacl::vector<double> vcl_eigenValues;
    viennacl::vector<double> vcl_basisCoeff;
    viennacl::matrix<double> vcl_velBasisField;
    viennacl::matrix<double> vcl_vortBasisField;
};

#endif
