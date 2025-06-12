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
    connect(ui->matrixValue, &QTextEdit::textChanged, this, &MainWindow::onMatrixValueChanged);
    connect(ui->zoomValue, &QLineEdit::textChanged, this, &MainWindow::onZoomValueChanged);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// Slot to update QTextEdit fields
void MainWindow::updateCameraValues(GLfloat* position) {
    QString matrixText;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j)
            matrixText += QString::number(position[i*4 + j], 'f', 4) + " ";
        matrixText += '\n';
    }

    ui->matrixValue->setText(matrixText.trimmed());
}

void MainWindow::updateCameraZoom(double zoom) {
    ui->zoomValue->setText(QString::number(zoom, 'f', 4));
}

void MainWindow::onZoomValueChanged() {
    QString value = ui->zoomValue->text();

    bool ok;
    double zoom = value.toDouble(&ok);

    if (ok) {
        ui->widget->setZoomLevel(zoom);
    } else {
        qDebug() << "Invalid zoom value!";
    }
}

void MainWindow::onMatrixValueChanged() {
    GLfloat newMatrix[16];
    if (setMatrixFromUI(newMatrix)) {
        ui->widget->updateViewMatrix(newMatrix);
    }
}

bool MainWindow::setMatrixFromUI(GLfloat* matrix) {
    QStringList values = ui->matrixValue->toPlainText().split(" ", Qt::SkipEmptyParts);

    if (values.size() != 16) {
        qDebug() << "Error: Matrix input must have exactly 16 values!";
        return false;
    }

    for (int i = 0; i < 16; ++i) {
        matrix[i] = values[i].toFloat();
    }

    return true;
}


