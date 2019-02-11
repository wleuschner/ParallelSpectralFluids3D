#include "Widgets/mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QGLFormat format;
    format = QGLFormat::defaultFormat();
    format.setProfile(QGLFormat::CoreProfile);
    format.setVersion(4,5);
    QGLFormat::setDefaultFormat(format);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
