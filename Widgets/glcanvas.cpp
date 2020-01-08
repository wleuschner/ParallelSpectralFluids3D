#include "glcanvas.h"
#include <iostream>
#include <QtOpenGL/QGLFramebufferObject>
#include <cstdlib>
#include "microprofile/microprofile.h"

#include <GL/gl.h>
#include <glm/gtc/matrix_transform.hpp>
#include <CL/cl.h>
#include <GL/glx.h>
#include <CL/cl_gl.h>

GLCanvas::GLCanvas(QWidget *parent) : QGLWidget(parent)
{
    meshVisible = true;
    voxelVisible = true;
    velocityVisible = false;
    particleVisible = true;
    psfSolver = NULL;
    solver = NULL;
    mesh = NULL;
    lifeTime = 30*60.0;
    imageNo = 0;
    record = false;
    benchmark = false;
    maxFrameNoBenchmark = 1000;
    frameNoBenchmark = 0;

    QGLFormat format = QGLFormat::defaultFormat();
    format.setProfile(QGLFormat::CoreProfile);
    format.setVersion(4,5);
    QGLFormat::setDefaultFormat(format);
    create();

    connect(&updateTimer,SIGNAL(timeout()),this,SLOT(simulate()));
    //connect(&updateTimer,SIGNAL(timeout()),this,SLOT(updateGL()));
    updateTimer.setInterval(1000/60);
    updateTimer.stop();
}

void GLCanvas::startBenchmark()
{
    benchmark = true;
    frameNoBenchmark = 0;
}

void GLCanvas::stopBenchmark()
{
    benchmark = false;
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
    updateGL();
}

void GLCanvas::showVoxel(bool state)
{
    voxelVisible = state;
    updateGL();
}

void GLCanvas::showVelocity(bool state)
{
    velocityVisible = state;
    updateGL();
}

void GLCanvas::showParticles(bool state)
{
    particleVisible = state;
    updateGL();
}

void GLCanvas::changeNumEigenfunctions(int n)
{
    updateTimer.stop();
    solver->setNumberEigenFunctions(n);
}

void GLCanvas::changeNumParticles(int n)
{
    updateTimer.stop();
    solver->changeNumParticles(n);
    updateTimer.start();
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

void GLCanvas::changeLifeTime(double val)
{
    solver->setLifeTime(val);
}

void GLCanvas::changeGravity(bool gravity)
{
    solver->setGravityActive(gravity);
}

void GLCanvas::changeGPU(bool gpu)
{
    if(gpu==true)
    {
        solver = psfSolverGPU;
    }
    else
    {
        solver = psfSolver;
    }
    solver->setMesh(mesh);
}

void GLCanvas::parameterChanged()
{

}

void GLCanvas::simulate()
{
    updateTimer.stop();
    //MicroProfileOnThreadCreate("Simulate");
    MICROPROFILE_SCOPEI("Canvas","simulate",MP_YELLOW);
    setUpdatesEnabled(false);
    makeCurrent();
    solver->integrate();
    updateGL();
    setUpdatesEnabled(true);
    doneCurrent();
    MicroProfileFlip(nullptr);
    updateTimer.start();
}

void GLCanvas::initializeGL()
{
    Display* display = glXGetCurrentDisplay();
    GLXContext gl_context = glXGetCurrentContext();

    // Get platform and device information
    cl_platform_id platform_id = NULL;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_GPU, 1,
            &device_id, &ret_num_devices);

    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM,(cl_context_properties) platform_id,
        CL_GLX_DISPLAY_KHR,(cl_context_properties) display,
        CL_GL_CONTEXT_KHR,(cl_context_properties) gl_context,
        0
    };

    cl_context_id = clCreateContext(props,1,&device_id,0,0,&ret);

    cl_queue = clCreateCommandQueue(cl_context_id,device_id,0,&ret);
    if(ret!=CL_SUCCESS) exit(-1);
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
    glDisable(GL_CULL_FACE);

    //Load Point Shader
    Shader pointVert(GL_VERTEX_SHADER,"Res/Effects/Points/points.vert");
    if(!pointVert.compile())
    {
        std::cout<<pointVert.compileLog()<<std::endl;
    }
    Shader pointFrag(GL_FRAGMENT_SHADER,"Res/Effects/Points/points.frag");
    if(!pointFrag.compile())
    {
        std::cout<<pointFrag.compileLog()<<std::endl;
    }
    pointsProgram = new ShaderProgram();
    pointsProgram->attachShader(pointVert);
    pointsProgram->attachShader(pointFrag);
    if(!pointsProgram->link())
    {
        std::cout<<pointsProgram->linkLog()<<std::endl;
    }
    pointsProgram->bind();

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

    //Load Volume Shader
    Shader volumeVert(GL_VERTEX_SHADER,"Res/Effects/Volume/volume.vert");
    if(!volumeVert.compile())
    {
        std::cout<<volumeVert.compileLog()<<std::endl;
    }
    Shader volumeFrag(GL_FRAGMENT_SHADER,"Res/Effects/Volume/volume.frag");
    if(!volumeFrag.compile())
    {
        std::cout<<volumeFrag.compileLog()<<std::endl;
    }
    volumeProgram = new ShaderProgram();
    volumeProgram->attachShader(volumeVert);
    volumeProgram->attachShader(volumeFrag);
    if(!volumeProgram->link())
    {
        std::cout<<volumeProgram->linkLog()<<std::endl;
    }
    volumeProgram->bind();

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

    glPointSize(2.0f);

    psfSolver = new PSFSolver();
    psfSolverGPU = new PSFSolverGPU(cl_context_id,device_id,cl_queue);
    psfSolverGPU->light = light;
    solver = psfSolverGPU;
}

void GLCanvas::paintGL()
{
    glViewport(0,0,width(),height());
    glClearColor(0.0,0.0,0.0,0.0);
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
            solver->drawGrid(lineProgram,pvm);
        }

        if(velocityVisible)
        {
            solver->drawVelocity(lineProgram,pvm);
        }

        if(particleVisible)
        {
            glDisable(GL_CULL_FACE);
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ONE,GL_ONE);

            Particle::setVertexAttribs();
            Particle::enableVertexAttribs();

            solver->view_mat = camera.getView();
            solver->camera_position = camera.getPosition();
            solver->drawParticles(volumeProgram,pvm);
            //solver->drawParticles(pointsProgram,pvm);
            glDisable(GL_BLEND);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_CULL_FACE);
        }

        if(meshVisible)
        {
            glDisable(GL_CULL_FACE);
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ONE,GL_ONE);


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
            Vertex::setVertexAttribs();
            Vertex::enableVertexAttribs();
            solver->drawParticles(lineProgram,pvm);
            glDisable(GL_BLEND);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_CULL_FACE);
        }

    }
    if(record)
    {
        grabFrameBuffer().save(QString("image")+QString::number(imageNo)+QString(".bmp"),"bmp");
        imageNo++;
    }
}

void GLCanvas::resizeGL(int w, int h)
{
    psfSolverGPU->viewport_size = glm::vec4(0.0,0.0,w,h);
    projection = glm::perspectiveFovRH(45.0f,(float)w,(float)h,0.1f,100.0f);
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
    simulate();
}
void GLCanvas::mouseMoveEvent(QMouseEvent *event)
{

}
