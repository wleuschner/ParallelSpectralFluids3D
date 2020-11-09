#include"abstractsolver.h"
#include<GL/glew.h>
#include<iostream>
#include<iomanip>
#include<chrono>

AbstractSolver::AbstractSolver()
{
    //maxParticles=1000;
    maxFramesBenchmark = 1000;
    isBenchmark = false;
    benchmarkFrameNo = 0;
    maxParticles=1000;
    gridVerts = NULL;
    gridIndices = NULL;
    velocityVerts = NULL;
    histogramTexture = NULL;
    volumeTextures[0] = NULL;
    volumeTextures[1] = NULL;

    for(unsigned int i=0;i<8;i++)
    {
        historyBuffer[i] = NULL;
        colorAttachments[i] = NULL;
        depthAttachments[i] = NULL;
    }
    mesh = NULL;
    gravityActive = false;
    resolution = 0.1;
    nEigenFunctions = 16;
    viscosity = 0.0f;
    timeStep = 1.0f/60.0f;
    lifeTime = 1.0*60.0;

    glm::vec4 verts[] = {glm::vec4(-1.0,-1.0,0.0,0.0),glm::vec4(-1.0,1.0,0.0,0.0),glm::vec4(1.0,1.0,0.0,0.0),glm::vec4(1.0,-1.0,0.0,0.0)};
    unsigned int indices[] = {0,1,2,2,0,3};
    glGenBuffers(1,&fullscreenVBO);
    glGenBuffers(1,&fullscreenIBO);
    glBindBuffer(GL_ARRAY_BUFFER,fullscreenVBO);
    glBufferData(GL_ARRAY_BUFFER,4*sizeof(glm::vec4),verts,GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,fullscreenIBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,6*sizeof(unsigned int),indices,GL_STATIC_DRAW);

    particles = new ParticleBuffer();
    particles->reserve(maxParticles);
}

void AbstractSolver::resize(unsigned int w,unsigned int h, bool gpu)
{
    if(gpu)
    {
        for(unsigned int i=0;i<8;i++)
        {
            if(historyBuffer[i]!=NULL)
            {
                delete historyBuffer[i];
            }
            if(colorAttachments[i]!=NULL)
            {
                delete colorAttachments[i];
            }
            if(depthAttachments[i]!=NULL)
            {
                delete depthAttachments[i];
            }
            colorAttachments[i] = new Texture();
            colorAttachments[i]->bind(0);
            colorAttachments[i]->createRenderImage(w,h);

            depthAttachments[i] = new Texture();
            depthAttachments[i]->bind(0);
            depthAttachments[i]->createDepthImage(w,h);

            historyBuffer[i] = new FrameBufferObject();
            historyBuffer[i]->bind();
            historyBuffer[i]->resize(w,h);
            historyBuffer[i]->attachColorImage(*colorAttachments[i],0);
            historyBuffer[i]->attachDepthImage(*depthAttachments[i]);
            if(!historyBuffer[i]->isComplete())
            {
                exit(-1);
            }
        }
    }
    else
    {
        glFramebufferTexture(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,0,0);
        glFramebufferTexture(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,0,0);
    }
}

void AbstractSolver::setMesh(Model* mesh)
{
    clearParticles();
    if(this->mesh!=NULL)
    {
        //delete this->mesh;
    }
    this->mesh = mesh;
    decMesh = mesh->voxelize(resolution);
    if(volumeTextures[0]!=NULL)
    {
        //volumeTextures[0]->destroy();
        delete volumeTextures[0];
    }
    if(volumeTextures[1]!=NULL)
    {
        //volumeTextures[1]->destroy();
        delete volumeTextures[1];
    }
    if(histogramTexture!=NULL)
    {
        //histogramTexture->destroy();
        delete histogramTexture;
    }
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

    glm::uvec3 dims = decMesh.getDimensions();

    volumeTextures[0] = new Texture3D();
    volumeTextures[0]->bind(0);
    volumeTextures[0]->createFloatRenderImage(256,256,256);
/*
    volumeTextures[1] = new Texture3D();
    volumeTextures[1]->bind(0);
    volumeTextures[1]->createFloatRenderImage(1024,1024,1024);*/

    histogramTexture = new Texture3D();
    histogramTexture->bind(0);
    histogramTexture->createRenderImage(256,256,256);

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

void AbstractSolver::setInitialVelocityField(const Eigen::VectorXf& field)
{
    basisCoeff = velBasisField.transpose()*field.cast<double>();
}

void AbstractSolver::setInitialVorticityField(const Eigen::VectorXf& field)
{
    basisCoeff = vortBasisField.transpose()*field.cast<double>();
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
    if(volumeTextures[0]!=NULL)
    {
        volumeTextures[0]->destroy();
        delete volumeTextures[0];
    }
    if(volumeTextures[1]!=NULL)
    {
        volumeTextures[1]->destroy();
        delete volumeTextures[1];
    }
    if(histogramTexture!=NULL)
    {
        histogramTexture->destroy();
        delete histogramTexture;
    }
    this->resolution = res;
    decMesh = mesh->voxelize(res);
    glm::uvec3 dims = decMesh.getDimensions();

    volumeTextures[0] = new Texture3D();
    volumeTextures[0]->bind(0);
    volumeTextures[0]->createFloatRenderImage(256,256,256);
/*
    volumeTextures[1] = new Texture3D();
    volumeTextures[1]->bind(0);
    volumeTextures[1]->createFloatRenderImage(512,512,512);*/

    histogramTexture = new Texture3D();
    histogramTexture->bind(0);
    histogramTexture->createRenderImage(256,256,256);

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

void AbstractSolver::changeNumParticles(unsigned int n)
{
    maxParticles = n;
    clearParticles();
}

void AbstractSolver::setLifeTime(float lt)
{
    lifeTime = lt*60.0;
}

void AbstractSolver::startBenchmark(bool benchmark)
{
    if(benchmark)
    {
        isBenchmark = true;
        benchmarkSums.resize(2);
        benchmarkSums[0] = 0.0;
        benchmarkSums[1] = 0.0;
        benchmarkFrameNo = 0;
        if(benchmark_file.is_open())
        {
            benchmark_file.close();
        }
        benchmark_file.open("benchmark.csv");
        benchmark_file.setf(std::ios::fixed,std::ios::floatfield);
        glm::uvec3 dims = decMesh.getDimensions();
        benchmark_file<<maxParticles<<std::endl;
        benchmark_file<<nEigenFunctions<<std::endl;
        benchmark_file<<dims.x<<","<<dims.y<<","<<dims.z<<std::endl;
        clearParticles();
        {
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            buildLaplace();
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            benchmark_file<<std::chrono::duration<double,std::milli>(end-begin).count()<<",";
        }
        {
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            buildEigenFunctions();
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            benchmark_file<<std::chrono::duration<double,std::milli>(end-begin).count()<<",";
        }
        {
            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            buildAdvection();
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            benchmark_file<<std::chrono::duration<double,std::milli>(end-begin).count()<<std::endl<<std::endl;
        }
    }
    else
    {
        benchmark_file<<std::endl<<std::endl;
        benchmark_file<<benchmarkSums[0]/maxFramesBenchmark<<","<<benchmarkSums[1]/maxFramesBenchmark;
        isBenchmark = false;
        benchmarkFrameNo = 0;
        if(benchmark_file.is_open())
        {
            benchmark_file.close();
        }
    }
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

const Eigen::VectorXf& AbstractSolver::getVelocityField()
{
    return velocityField;
}

const Eigen::VectorXf& AbstractSolver::getVorticityField()
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

void AbstractSolver::beginBenchmark()
{

}

void AbstractSolver::endBenchmark()
{

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
