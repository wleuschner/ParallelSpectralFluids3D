#ifndef __PSF_SOLVER_GPU_H_
#define __PSF_SOLVER_GPU_H_
#include "../abstractsolver.h"
#include <viennacl/vector.hpp>
#include <viennacl/matrix.hpp>
#include <viennacl/compressed_matrix.hpp>

class PSFSolverGPU : public AbstractSolver
{
public:
    PSFSolverGPU();
    void integrate();
protected:
    void buildLaplace();
    void buildAdvection();
private:
    viennacl::compressed_matrix<double> vcl_curl;

    std::vector<viennacl::vector<double>> vcl_eigenFunctions;

    viennacl::vector<double> vcl_vorticityField;
    viennacl::vector<double> vcl_velocityField;

    std::vector<viennacl::matrix<double>> vcl_advection;

    viennacl::vector<double> vcl_gravity;
    viennacl::vector<double> vcl_eigenValues;
    viennacl::vector<double> vcl_basisCoeff;
    viennacl::matrix<double> vcl_velBasisField;
    viennacl::matrix<double> vcl_vortBasisField;

    std::vector<Particle> particles;
};

#endif
