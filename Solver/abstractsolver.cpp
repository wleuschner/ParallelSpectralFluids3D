#include"abstractsolver.h"
#include<GL/gl.h>
#include<iostream>

AbstractSolver::AbstractSolver()
{
    gridVerts = NULL;
    gridIndices = NULL;
    velocityVerts = NULL;
    mesh = NULL;
    resolution = 0.2;
    nEigenFunctions = 16;
    viscosity = 0.0f;
    timeStep = 1.0f/60.0f;
}

void AbstractSolver::setMesh(Model* mesh)
{
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
    particleVerts = new VertexBuffer();
    std::vector<unsigned int> indices(decMesh.getNumEdges());
    unsigned int e=0;
    for(EdgeIterator it=decMesh.getEdgeIteratorBegin();it!=decMesh.getEdgeIteratorEnd();++it)
    {
        if(it->inside==GridState::INSIDE)
        {
            indices[e] = it->v1;
            e++;
            indices[e] = it->v2;
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
    this->nEigenFunctions = n;
    buildLaplace();
    buildEigenFunctions();
    buildAdvection();
}

void AbstractSolver::setResolution(double res)
{
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

const std::vector<glm::dvec3>& AbstractSolver::getParticles()
{
    return particles;
}

void AbstractSolver::clearParticles()
{
    particles.clear();
}

void AbstractSolver::addParticle(glm::dvec3 particle)
{
    if(decMesh.isPointInside(particle))
    {
        particles.push_back(particle);
    }
}

unsigned int AbstractSolver::getNumParticles()
{
    return particles.size();
}

void AbstractSolver::drawGrid(ShaderProgram* program,const glm::mat4& pvm)
{
    gridIndices->bind();
    gridVerts->bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
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



            double s1 = decMesh.getFaceSignum(it->f1,it->v1,it->v2,it->v6,it->v5);
            double s2 = decMesh.getFaceSignum(it->f2,it->v3,it->v4,it->v8,it->v7);
            double s3 = decMesh.getFaceSignum(it->f3,it->v7,it->v8,it->v5,it->v6);
            double s4 = decMesh.getFaceSignum(it->f4,it->v4,it->v3,it->v2,it->v1);
            double s5 = decMesh.getFaceSignum(it->f5,it->v4,it->v1,it->v5,it->v8);
            double s6 = decMesh.getFaceSignum(it->f6,it->v2,it->v3,it->v7,it->v6);

            glm::vec3 center = vertices[it->v5].pos+0.5f*(vertices[it->v3].pos-vertices[it->v5].pos);
            glm::vec3 velDir = 0.5f*glm::vec3((velocityField(f1.id)*f1.normal+velocityField(f2.id)*f2.normal))+
                               0.5f*glm::vec3((velocityField(f3.id)*f3.normal+velocityField(f4.id)*f4.normal))+
                               0.5f*glm::vec3((velocityField(f5.id)*f5.normal+velocityField(f6.id)*f6.normal));

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
    glDrawArrays(GL_LINES,0,decMesh.getNumVoxels()*2);
}

void AbstractSolver::drawParticles(ShaderProgram* program,const glm::mat4& pvm)
{
    particleVerts->bind();
    std::vector<Vertex> particleV(getNumParticles());
    for(unsigned int i=0;i<getNumParticles();i++)
    {
        particleV[i].pos = glm::vec3(particles[i]);
    }
    particleVerts->upload(particleV);
    particleVerts->bind();
    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();
    program->bind();
    program->uploadMat4("pvm",pvm);
    glDrawArrays(GL_POINTS,0,getNumParticles());
}

void AbstractSolver::buildEigenFunctions()
{

}
