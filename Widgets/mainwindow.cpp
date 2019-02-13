#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <iostream>
#include <QPainter>
#include "../DEC/dec.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    connect(ui->actionOpen,SIGNAL(triggered()),this,SLOT(openFile()));
    connect(ui->chkShowVoxel,SIGNAL(toggled(bool)),ui->canvas,SLOT(showVoxel(bool)));
    connect(ui->chkShowMesh,SIGNAL(toggled(bool)),ui->canvas,SLOT(showMesh(bool)));
    connect(ui->chkVelocity,SIGNAL(toggled(bool)),ui->canvas,SLOT(showVelocity(bool)));
    connect(ui->spinNumEigenFunctions,SIGNAL(valueChanged(int)),ui->canvas,SLOT(changeNumEigenfunctions(int)));
    connect(ui->spinVoxelSize,SIGNAL(valueChanged(double)),ui->canvas,SLOT(changeResolution(double)));
    connect(ui->spinViscosity,SIGNAL(valueChanged(double)),ui->canvas,SLOT(changeViscosity(double)));
    connect(ui->spinTimestep,SIGNAL(valueChanged(double)),ui->canvas,SLOT(changeTimestep(double)));
    connect(ui->spinLifeTime,SIGNAL(valueChanged(double)),ui->canvas,SLOT(changeLifeTime(double)));

    mesh = NULL;
}

MainWindow::~MainWindow()
{
    delete ui;
    if(mesh!=NULL)
    {
        delete mesh;
    }
}

void MainWindow::openFile()
{
    QString fileName = QFileDialog::getOpenFileName(this,"Mesh File");
    if(!fileName.isEmpty())
    {
        //Init Mesh
        if(mesh!=NULL)
        {
            mesh->release();
            delete mesh;
        }
        mesh = new Model();
        if(!mesh->load(fileName.toStdString()))
        {
            delete mesh;
            mesh = NULL;
        }

        ui->canvas->setMesh(mesh);
        //ui->renderWidget->setMinimumSize(originalImage.width(),originalImage.height());
        //ui->renderWidget->setPixmap(QPixmap::fromImage(originalImage));

    }
}
