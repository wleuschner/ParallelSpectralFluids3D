#include"abstractsolver.h"
#include<GL/gl.h>
#include<iostream>

AbstractSolver::AbstractSolver()
{
    //maxParticles=1000;
    maxParticles=400000;
    gridVerts = NULL;
    gridIndices = NULL;
    velocityVerts = NULL;
    mesh = NULL;
    gravityActive = false;
    resolution = 0.2;
    nEigenFunctions = 16;
    viscosity = 0.0f;
    timeStep = 1.0f/60.0f;

    particles = new ParticleBuffer();
    particles->reserve(maxParticles);
}

void AbstractSolver::setMesh(Model* mesh)
{
    clearParticles();
    if(this->mesh!=NULL)
    {
        delete this->mesh;
    }
    this->mesh = mesh;
    decMesh = mesh->voxelize(resolution);

    if(gridVerts!=NULL)
    {
        delete gridVerts;
    }
    if(gridIndices!=NULL)
    {
        delete gridIndices;
    }
    if(velocityVerts!=NULL)
    {
        delete velocityVerts;
    }
    gridVerts = new VertexBuffer();
    gridIndices = new IndexBuffer();
    velocityVerts = new VertexBuffer();
    std::vector<unsigned int> indices(2*decMesh.getNumEdges());
    unsigned int e=0;
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();++it)
    {
        if(it->inside==GridState::INSIDE)
        {
            indices[e] = decMesh.signedIdToIndex(it->v1);
            e++;
            indices[e] = decMesh.signedIdToIndex(it->v2);
            e++;
        }
    }
    gridIndices->bind();
    gridIndices->upload(indices);
    gridVerts->bind();
    gridVerts->upload(mesh->getVertices());

    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setInitialVelocityField(const Eigen::VectorXd& field)
{
    basisCoeff = velBasisField.transpose()*field;
}

void AbstractSolver::setInitialVorticityField(const Eigen::VectorXd& field)
{
    basisCoeff = vortBasisField.transpose()*field;
}

void AbstractSolver::setNumberEigenFunctions(unsigned int n)
{
    clearParticles();
    this->nEigenFunctions = n;
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setResolution(double res)
{
    clearParticles();
    this->resolution = res;
    decMesh = mesh->voxelize(res);
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setViscosity(double visc)
{
    this->viscosity = visc;
}

void AbstractSolver::setTimestep(double timestep)
{
    this->timeStep = timestep;
}

void AbstractSolver::setGravityActive(bool state)
{
    gravityActive = state;
}

Model* AbstractSolver::getMesh()
{
    return mesh;
}

DECMesh3D& AbstractSolver::getDECMesh()
{
    return decMesh;
}

unsigned int AbstractSolver::getNumEigenFunctions()
{
    return nEigenFunctions;
}

double AbstractSolver::getTimestep()
{
    return timeStep;
}

const Eigen::VectorXd& AbstractSolver::getEigenFunction(unsigned int n)
{
    return eigenFunctions[n];
}

const Eigen::MatrixXd& AbstractSolver::getVelocityBasisField()
{
    return velBasisField;
}

const Eigen::MatrixXd& AbstractSolver::getVorticityBasisField()
{
    return vortBasisField;
}

const Eigen::VectorXd& AbstractSolver::getBasisCoefficients()
{
    return basisCoeff;
}

const Eigen::VectorXd& AbstractSolver::getVelocityField()
{
    return velocityField;
}

const Eigen::VectorXd& AbstractSolver::getVorticityField()
{
    return vorticityField;
}

double AbstractSolver::getMaxVorticity()
{
    return maxRotation;
}

double AbstractSolver::getMinVorticity()
{
    return minRotation;
}

const std::vector<Particle>& AbstractSolver::getParticles()
{
    return particles->getParticles();
}

void AbstractSolver::clearParticles()
{
    simTime = 0.0;
    particlePointer = 0;
    particles->clear();
    particles->reserve(maxParticles);
}

void AbstractSolver::addParticle(Particle particle)
{
    std::vector<Particle>& parts = particles->getParticles();
    if(decMesh.isPointInside(particle.position)&&
       simTime>=parts[particlePointer].position.w)
    {
        parts[particlePointer].position = particle.position;
        parts[particlePointer].position.w = simTime+particle.position.w;
        particlePointer = (particlePointer+1)%maxParticles;
    }
}

unsigned int AbstractSolver::getNumParticles()
{
    return particles->getNumParticles();
}

void AbstractSolver::drawGrid(ShaderProgram* program,const glm::mat4& pvm)
{
    gridIndices->bind();
    gridVerts->bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    program->uploadVec4("color",glm::vec4(1.0,0.0,0.0,1.0));
    glDrawElements(GL_LINES,decMesh.getNumEdges()*2,GL_UNSIGNED_INT,(void*)0);
}

void AbstractSolver::drawVelocity(ShaderProgram* program,const glm::mat4& pvm)
{
    velocityVerts->bind();
    unsigned int e=0;
    std::vector<Vertex> vertices = mesh->getVertices();
    std::vector<Vertex> velVerts(decMesh.getNumVoxels()*2);
    for(VoxelIterator it=decMesh.getVoxelIteratorBegin();it<decMesh.getVoxelIteratorEnd();++it)
    {
        if(it->inside==GridState::INSIDE)
        {
            Face3D f1 = decMesh.getFace(it->f1);
            Face3D f2 = decMesh.getFace(it->f2);
            Face3D f3 = decMesh.getFace(it->f3);
            Face3D f4 = decMesh.getFace(it->f4);
            Face3D f5 = decMesh.getFace(it->f5);
            Face3D f6 = decMesh.getFace(it->f6);

            double s1 = decMesh.getFaceSignum(it->f1);
            double s2 = decMesh.getFaceSignum(it->f2);
            double s3 = decMesh.getFaceSignum(it->f3);
            double s4 = decMesh.getFaceSignum(it->f4);
            double s5 = decMesh.getFaceSignum(it->f5);
            double s6 = decMesh.getFaceSignum(it->f6);

            glm::vec3 center = (f1.center+f2.center)/2.0;
            glm::vec3 velDir = 0.5f*glm::vec3((s5*velocityField(f5.id)*f5.normal-s6*velocityField(f6.id)*f6.normal))+
                               0.5f*glm::vec3((s3*velocityField(f3.id)*f3.normal-s4*velocityField(f4.id)*f4.normal))+
                               0.5f*glm::vec3((s1*velocityField(f1.id)*f1.normal-s2*velocityField(f2.id)*f2.normal));

            velVerts[e].pos = center;
            e++;
            velVerts[e].pos = center+velDir;
            e++;
        }
    }
    velocityVerts->upload(velVerts);
    velocityVerts->bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    program->uploadVec4("color",glm::vec4(1.0,0.0,0.0,1.0));
    glDrawArrays(GL_LINES,0,decMesh.getNumVoxels()*2);
}

void AbstractSolver::buildEigenFunctions()
{

}
