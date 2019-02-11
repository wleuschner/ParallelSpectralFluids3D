#ifndef __PSF_SOLVER_GPU_H_
#define __PSF_SOLVER_GPU_H_
#include"../abstractsolver.h"

class PSFSolverGPU : public AbstractSolver
{
public:
    PSFSolverGPU();
    void integrate();
};

#endif
