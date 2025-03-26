#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QtDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // Connect GLWidget's signal to update QLineEdit fields
    connect(ui->widget, &GLWidget::cameraUpdated, this, &MainWindow::updateCameraValues);
    connect(ui->widget, &GLWidget::cameraUpdatedZoom, this, &MainWindow::updateCameraZoom);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Slot to update QLineEdit fields
void MainWindow::updateCameraValues(const QVector3D &position) {
    ui->xValue->setText(QString::number(position.x(), 'f', 4));
    ui->yValue->setText(QString::number(position.y(), 'f', 4));
    ui->zValue->setText(QString::number(position.z(), 'f', 4));
}

void MainWindow::updateCameraZoom(double zoom) {
    ui->zoomValue->setText(QString::number(zoom, 'f', 4));
}
