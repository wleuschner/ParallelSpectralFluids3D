#ifndef GLCanvas_H
#define GLCanvas_H

#include <GL/glew.h>
#include <QObject>
#include <QtOpenGL/QGLWidget>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QTimer>
#include "../Graphics/Camera/Camera.h"
#include "../Graphics/Model/Model.h"
#include "../Solver/abstractsolver.h"
#include "../Solver/SolverImpl/psfsolver.h"
#include "../Solver/SolverImpl/psfsolvergpu.h"



class GLCanvas : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLCanvas(QWidget *parent = nullptr);
    void setMesh(Model* mesh);
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);
    void keyPressEvent(QKeyEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
signals:

public slots:
    void startBenchmark();
    void stopBenchmark();
    void showMesh(bool state);
    void showVoxel(bool state);
    void showVelocity(bool state);
    void showParticles(bool state);
    void changeNumEigenfunctions(int n);
    void changeResolution(double resolution);
    void changeViscosity(double visc);
    void changeTimestep(double val);
    void changeLifeTime(double val);
    void changeNumParticles(int val);
    void changeGravity(bool gravity);
    void changeGPU(bool gpu);
    void parameterChanged();
    void simulate();
private:
    QTimer updateTimer;
    cl_context cl_context_id;
    cl_command_queue cl_queue;
    cl_device_id device_id;

    bool benchmark;
    unsigned int frameNoBenchmark;
    unsigned int maxFrameNoBenchmark;

    bool record;
    unsigned int imageNo;
    double lifeTime;

    PSFSolver* psfSolver;
    PSFSolverGPU* psfSolverGPU;
    AbstractSolver* solver;
    bool meshVisible;
    bool voxelVisible;
    bool velocityVisible;
    bool particleVisible;
    Camera camera;
    Model* mesh;
    unsigned int vao;
    glm::mat4 projection;
    ShaderProgram* volumeProgram;
    ShaderProgram* phongProgram;
    ShaderProgram* lineProgram;
    ShaderProgram* pointsProgram;
    Light light;
};

#endif // GLCanvas_H
