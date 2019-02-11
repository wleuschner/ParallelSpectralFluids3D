#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <Eigen/Eigen>
#include <assimp/cimport.h>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "Widgets/glcanvas.h"
#include "Graphics/Model/Model.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
private slots:
    void openFile();

protected:
private:
    void voxelize(unsigned int w,unsigned int h);
    bool checkVoxel(unsigned char* offs,int w,int h);

    Model* mesh;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
