#ifndef __ABSTRACT_SOLVER_H_
#define __ABSTRACT_SOLVER_H_
#include <Eigen/Eigen>
#include "../Graphics/Model/Model.h"
#include "../Graphics/ParticleBuffer/ParticleBuffer.h"
#include "../Graphics/Particle/Particle.h"
#include "../DEC/decmesh3d.h"

class AbstractSolver
{
public:
    AbstractSolver();

    virtual void integrate()=0;
    void setMesh(Model* mesh);

    void setInitialVelocityField(const Eigen::VectorXd& field);
    void setInitialVorticityField(const Eigen::VectorXd& field);

    void setNumberEigenFunctions(unsigned int n);
    void setResolution(double res);
    void setViscosity(double visc);
    void setTimestep(double timestep);
    void setGravityActive(bool state);

    Model* getMesh();
    DECMesh3D& getDECMesh();

    unsigned int getNumEigenFunctions();
    double getTimestep();
    const Eigen::VectorXd& getEigenFunction(unsigned int n);
    const Eigen::MatrixXd& getVelocityBasisField();
    const Eigen::MatrixXd& getVorticityBasisField();
    const Eigen::VectorXd& getBasisCoefficients();
    const Eigen::VectorXd& getVelocityField();
    const Eigen::VectorXd& getVorticityField();
    double getMaxVorticity();
    double getMinVorticity();

    const std::vector<Particle>& getParticles();
    void clearParticles();
    void addParticle(Particle particle);
    unsigned int getNumParticles();

    void drawGrid(ShaderProgram* program,const glm::mat4& pvm);
    void drawVelocity(ShaderProgram* program,const glm::mat4& pvm);
    virtual void drawParticles(ShaderProgram* program,const glm::mat4& pvm) = 0;

protected:
    unsigned int maxParticles;
    unsigned int particlePointer;
    double simTime;

    virtual void buildLaplace()=0;
    virtual void buildAdvection()=0;

    void buildEigenFunctions();

    VertexBuffer* gridVerts;
    IndexBuffer* gridIndices;

    VertexBuffer* velocityVerts;

    double minRotation;
    double maxRotation;

    unsigned int nEigenFunctions;
    double resolution;
    double timeStep;
    double viscosity;

    DECMesh3D decMesh;
    Model* mesh;

    bool gravityActive;

    Eigen::SparseMatrix<double> curl;

    std::vector<Eigen::VectorXd> eigenFunctions;

    Eigen::VectorXd vorticityField;
    Eigen::VectorXd velocityField;

    std::vector<Eigen::MatrixXd> advection;

    Eigen::VectorXd gravity;
    Eigen::VectorXd eigenValues;
    Eigen::VectorXd basisCoeff;
    Eigen::MatrixXd velBasisField;
    Eigen::MatrixXd vortBasisField;

    ParticleBuffer* particles;
};

#endif
