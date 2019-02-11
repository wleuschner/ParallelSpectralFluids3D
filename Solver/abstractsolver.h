#ifndef __ABSTRACT_SOLVER_H_
#define __ABSTRACT_SOLVER_H_
#include <Eigen/Eigen>
#include "../Graphics/Model/Model.h"
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

    const std::vector<glm::dvec3>& getParticles();
    void clearParticles();
    void addParticle(glm::dvec3 particle);
    unsigned int getNumParticles();

    void drawGrid(ShaderProgram* program,const glm::mat4& pvm);
    void drawVelocity(ShaderProgram* program,const glm::mat4& pvm);
    void drawParticles(ShaderProgram* program,const glm::mat4& pvm);

protected:
    virtual void buildLaplace()=0;
    virtual void buildAdvection()=0;

    void buildEigenFunctions();

    VertexBuffer* gridVerts;
    IndexBuffer* gridIndices;

    VertexBuffer* velocityVerts;
    VertexBuffer* particleVerts;

    double minRotation;
    double maxRotation;

    unsigned int nEigenFunctions;
    double resolution;
    double timeStep;
    double viscosity;

    DECMesh3D decMesh;
    Model* mesh;

    Eigen::SparseMatrix<double> curl;

    std::vector<Eigen::VectorXd> eigenFunctions;

    Eigen::VectorXd vorticityField;
    Eigen::VectorXd velocityField;

    std::vector<Eigen::MatrixXd> advection;

    Eigen::VectorXd eigenValues;
    Eigen::VectorXd basisCoeff;
    Eigen::MatrixXd velBasisField;
    Eigen::MatrixXd vortBasisField;

    std::vector<glm::dvec3> particles;
};

#endif
