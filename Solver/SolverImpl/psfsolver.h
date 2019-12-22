#ifndef __PSF_SOLVER_H_
#define __PSF_SOLVER_H_
#include"../abstractsolver.h"
#include"../../Graphics/Vertex/Vertex.h"
#include"../../Graphics/VertexBuffer/VertexBuffer.h"
#include"../../Graphics/IndexBuffer/IndexBuffer.h"

class PSFSolver : public AbstractSolver
{
public:
    PSFSolver();
    void integrate();

    void drawParticles(ShaderProgram* program,const glm::mat4& pvm);

protected:
    void buildLaplace();
    void buildAdvection();
};

#endif
