#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGLWidget>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void onMatrixValueChanged();
    void onZoomValueChanged();

private slots:
    void updateCameraValues(GLfloat* position);
    void updateCameraZoom(double zoom);

private:
    Ui::MainWindow *ui;
    bool setMatrixFromUI(GLfloat* matrix);
};

#endif // MAINWINDOW_H
