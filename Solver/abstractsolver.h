#ifndef __ABSTRACT_SOLVER_H_
#define __ABSTRACT_SOLVER_H_
#include <Eigen/Eigen>
#include <fstream>
#include "../Graphics/Model/Model.h"
#include "../Graphics/ParticleBuffer/ParticleBuffer.h"
#include "../Graphics/Particle/Particle.h"
#include "../Graphics/Texture3D/Texture3D.h"
#include "../Graphics/SSBO/SSBO.h"
#include "../Graphics/FrameBufferObject/FrameBufferObject.h"
#include "../DEC/decmesh3d.h"

class AbstractSolver
{
public:
    AbstractSolver();

    virtual void integrate()=0;
    void setMesh(Model* mesh);

    void setInitialVelocityField(const Eigen::VectorXf& field);
    void setInitialVorticityField(const Eigen::VectorXf& field);

    void resize(unsigned int w,unsigned int h);
    void setNumberEigenFunctions(unsigned int n);
    void changeNumParticles(unsigned int n);
    void setResolution(double res);
    void setViscosity(double visc);
    void setTimestep(double timestep);
    void setGravityActive(bool state);
    void setLifeTime(float lt);
    void startBenchmark(bool benchmark);

    void beginBenchmark();
    void endBenchmark();

    Model* getMesh();
    DECMesh3D& getDECMesh();

    unsigned int getNumEigenFunctions();
    double getTimestep();
    const Eigen::VectorXd& getEigenFunction(unsigned int n);
    const Eigen::MatrixXd& getVelocityBasisField();
    const Eigen::MatrixXd& getVorticityBasisField();
    const Eigen::VectorXd& getBasisCoefficients();
    const Eigen::VectorXf& getVelocityField();
    const Eigen::VectorXf& getVorticityField();
    double getMaxVorticity();
    double getMinVorticity();

    const std::vector<Particle>& getParticles();
    void clearParticles();
    void addParticle(Particle particle);
    unsigned int getNumParticles();

    void drawGrid(ShaderProgram* program,const glm::mat4& pvm);
    void drawVelocity(ShaderProgram* program,const glm::mat4& pvm);
    virtual void drawParticles(ShaderProgram* program,const glm::mat4& pvm) = 0;

    glm::vec4 viewport_size;
    glm::vec3 camera_position;
    glm::mat4 view_mat;

protected:
    unsigned int maxFramesBenchmark;
    bool isBenchmark;
    unsigned int benchmarkFrameNo;
    std::ofstream benchmark_file;
    std::vector<double> benchmarkSums;

    unsigned int maxParticles;
    unsigned int particlePointer;
    double simTime;

    virtual void buildLaplace()=0;
    virtual void buildAdvection()=0;

    void buildEigenFunctions();

    VertexBuffer* gridVerts;
    IndexBuffer* gridIndices;

    Texture3D* histogramTexture;
    Texture3D* volumeTextures[2];
    VertexBuffer* velocityVerts;

    unsigned int currentHistory=0;
    FrameBufferObject* historyBuffer[8];
    Texture* colorAttachments[8];
    Texture* depthAttachments[8];

    double minRotation;
    double maxRotation;

    float lifeTime;
    unsigned int nEigenFunctions;
    double resolution;
    double timeStep;
    double viscosity;

    DECMesh3D decMesh;
    Model* mesh;

    bool gravityActive;

    unsigned int fullscreenVBO;
    unsigned int fullscreenIBO;

    Eigen::SparseMatrix<double> curl;

    std::vector<Eigen::VectorXd> eigenFunctions;

    Eigen::VectorXf vorticityField;
    Eigen::VectorXf velocityField;

    std::vector<Eigen::MatrixXd> advection;

    Eigen::VectorXd gravity;
    Eigen::VectorXd eigenValues;
    Eigen::VectorXd basisCoeff;
    Eigen::MatrixXd velBasisField;
    Eigen::MatrixXd vortBasisField;

    ParticleBuffer* particles;
};

#endif
