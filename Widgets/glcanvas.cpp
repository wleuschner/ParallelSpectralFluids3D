#include "glcanvas.h"
#include <iostream>
#include <QtOpenGL/QGLFramebufferObject>
#include <cstdlib>

#include <GL/gl.h>
#include <glm/gtc/matrix_transform.hpp>

GLCanvas::GLCanvas(QWidget *parent) : QGLWidget(parent)
{
    meshVisible = true;
    voxelVisible = true;
    velocityVisible = false;
    psfSolver = NULL;
    solver = NULL;
    mesh = NULL;
    imageNo = 0;
    record = false;
    QGLFormat format = QGLFormat::defaultFormat();
    format.setProfile(QGLFormat::CoreProfile);
    format.setVersion(4,5);
    QGLFormat::setDefaultFormat(format);
    create();

    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(simulate()));
    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(updateGL()));
    updateTimer.setInterval(1000/60);
    updateTimer.stop();
}

void GLCanvas::setMesh(Model* mesh)
{
    updateTimer.stop();
    makeCurrent();
    this->mesh = mesh;
    solver->setMesh(mesh);
    QGLFramebufferObject::bindDefault();
    doneCurrent();
}

void GLCanvas::showMesh(bool state)
{
    meshVisible = state;
    update();
}

void GLCanvas::showVoxel(bool state)
{
    voxelVisible = state;
    update();
}

void GLCanvas::changeNumEigenfunctions(int n)
{
    updateTimer.stop();
    solver->setNumberEigenFunctions(n);
}

void GLCanvas::changeResolution(double resolution)
{
    updateTimer.stop();
    solver->setResolution(resolution);
}

void GLCanvas::changeViscosity(double visc)
{
    solver->setViscosity(visc);
}

void GLCanvas::changeTimestep(double val)
{
    solver->setTimestep(val);
}


void GLCanvas::parameterChanged()
{

}

void GLCanvas::simulate()
{
    if(solver->getNumParticles()<1000000)
    {
        solver->addParticle(glm::dvec3((rand()%1024)/1024.0-0.5,-0.9,(rand()%1024)/1024.0-0.5));
    }
    solver->integrate();
}

void GLCanvas::initializeGL()
{
    glewExperimental = true;
    if(glewInit()!=GLEW_OK)
    {
        std::cout<<"GLEW ERROR"<<std::endl;
    }
    glGenVertexArrays(1,&vao);
    glBindVertexArray(vao);

    glClearColor(1.0,0.0,0.0,1.0);
    glClearDepth(1.0f);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_VERTEX_ARRAY);
    glDisable(GL_CULL_FACE);
    //Load Line Shader
    Shader lineVert(GL_VERTEX_SHADER,"Res/Effects/Line/line.vert");
    if(!lineVert.compile())
    {
        std::cout<<lineVert.compileLog()<<std::endl;
    }
    Shader lineFrag(GL_FRAGMENT_SHADER,"Res/Effects/Line/line.frag");
    if(!lineFrag.compile())
    {
        std::cout<<lineFrag.compileLog()<<std::endl;
    }
    lineProgram = new ShaderProgram();
    lineProgram->attachShader(lineVert);
    lineProgram->attachShader(lineFrag);
    if(!lineProgram->link())
    {
        std::cout<<lineProgram->linkLog()<<std::endl;
    }
    lineProgram->bind();

    //Load Phong Shader
    Shader phongVert(GL_VERTEX_SHADER,"Res/Effects/Phong/phong.vert");
    if(!phongVert.compile())
    {
        std::cout<<phongVert.compileLog()<<std::endl;
    }
    Shader phongFrag(GL_FRAGMENT_SHADER,"Res/Effects/Phong/phong.frag");
    if(!phongFrag.compile())
    {
        std::cout<<phongFrag.compileLog()<<std::endl;
    }
    phongProgram = new ShaderProgram();
    phongProgram->attachShader(phongVert);
    phongProgram->attachShader(phongFrag);
    if(!phongProgram->link())
    {
        std::cout<<phongProgram->linkLog()<<std::endl;
    }
    phongProgram->bind();

    light = Light(glm::vec3(20.0,20.0,0.0));

    Vertex::setVertexAttribs();
    Vertex::enableVertexAttribs();

    psfSolver = new PSFSolver();
    solver = psfSolver;

}

void GLCanvas::paintGL()
{
    glViewport(0,0,width(),height());
    glClearColor(1.0,1.0,1.0,1.0);
    //glClearColor(0.0,0.0,0.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    if(mesh!=NULL)
    {
        glm::mat4 view = camera.getView();
        glm::mat4 model = mesh->getModelMat();
        glm::mat4 modelView = view*mesh->getModelMat();
        glm::mat4 pvm = projection*modelView;

        if(voxelVisible)
        {
            psfSolver->drawGrid(lineProgram,pvm);
        }

        if(velocityVisible)
        {
            psfSolver->drawVelocity(lineProgram,pvm);
        }

        if(meshVisible)
        {
            glDisable(GL_CULL_FACE);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ZERO,GL_SRC_COLOR);


            mesh->bind();
            Vertex::setVertexAttribs();
            Vertex::enableVertexAttribs();
            glm::mat3 normalMatrix = glm::mat3(glm::transpose(glm::inverse((modelView))));
            phongProgram->bind();
            phongProgram->uploadMat4("modelview",modelView);
            phongProgram->uploadMat4("pvm",pvm);
            phongProgram->uploadMat4("view",view);
            phongProgram->uploadMat3("normalMatrix",normalMatrix);
            //phongProgram->uploadVec3("cPos",camera.getPosition());
            phongProgram->uploadLight("light0",light,view);
            glDrawElements(GL_TRIANGLES,mesh->getIndices().size(),GL_UNSIGNED_INT,(void*)0);
            glDisable(GL_BLEND);
            glEnable(GL_CULL_FACE);
        }

        psfSolver->drawParticles(lineProgram,pvm);
    }
    if(record)
    {
        grabFrameBuffer().save(QString("image")+QString::number(imageNo)+QString(".bmp"),"bmp");
        imageNo++;
    }
}

void GLCanvas::resizeGL(int w, int h)
{
    projection = glm::perspectiveFov(45.0f,(float)w,(float)h,0.1f,100.0f);
}

void GLCanvas::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_W:
        camera.translate(0.1f*camera.getForwardVec());
        break;
    case Qt::Key_S:
        camera.translate(-0.1f*camera.getForwardVec());
        break;
    case Qt::Key_A:
        camera.translate(-0.1f*camera.getStrafeVec());
        break;
    case Qt::Key_D:
        camera.translate(0.1f*camera.getStrafeVec());
        break;
    case Qt::Key_PageUp:
        camera.translate(0.1f*camera.getUpVector());
        break;
    case Qt::Key_PageDown:
        camera.translate(-0.1f*camera.getUpVector());
        break;
    case Qt::Key_Up:
        camera.rotate(0.1f,camera.getStrafeVec());
        break;
    case Qt::Key_Down:
        camera.rotate(-0.1f,camera.getStrafeVec());
        break;
    case Qt::Key_Left:
        camera.rotate(0.1f,camera.getUpVector());
        break;
    case Qt::Key_Right:
        camera.rotate(-0.1f,camera.getUpVector());
        break;
    case Qt::Key_Plus:
        simulate();
        updateGL();
        break;
    case Qt::Key_Space:
        if(updateTimer.isActive())
        {
            updateTimer.stop();
        }
        else
        {
            updateTimer.start();
        }
        break;
    case Qt::Key_V:
        if(record)
        {
            record = false;
            imageNo = 0;
        }
        else
        {
            record = true;
        }
        break;
    }
    update();
}
void GLCanvas::mouseMoveEvent(QMouseEvent *event)
{

}
