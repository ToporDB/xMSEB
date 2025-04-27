#include <QtWidgets/QApplication>
#include <QFile>
#include <QTextStream>
#include <QLineEdit>
#include <QTextEdit>
#include <QDebug>
#include "mainwindow.h"


QString arg(QString argname) {
    int pos = qApp->arguments().indexOf(QRegExp("-" + argname));
    if (pos != -1 && pos + 1 < qApp->arguments().size()) {
        return qApp->arguments().at(pos + 1);
    }
    return "";
}

bool loadViewMatrixAndZoom(const QString& filename, double& zoom, QString& matrixString) {
    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qWarning() << "Failed to open view matrix file:" << filename;
        return false;
    }

    QTextStream in(&file);
    QVector<float> values;
    while (!in.atEnd()) {
        QString line = in.readLine();
        for (const QString& numStr : line.split(QRegExp("\\s+"), Qt::SkipEmptyParts)) {
            values.append(numStr.toFloat());
        }
    }
    file.close();

    if (values.size() != 17) {
        qWarning() << "Invalid view matrix file: Expected 17 floats, got" << values.size();
        return false;
    }

    zoom = static_cast<double>(values[0]);

    // Build the matrix string for QTextEdit
    QStringList matrixParts;
    for (int i = 1; i <= 16; ++i) {
        matrixParts << QString::number(values[i], 'f', 6);
    }
    matrixString = matrixParts.join(" ");

    return true;
}

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow w;
    w.setWindowTitle("Fibviewer");

    // Check if -viewmatrix is passed
    QString viewMatrixFile = arg("viewmatrix");
    if (!viewMatrixFile.isEmpty()) {
        double zoom = 1.0;
        QString matrixString;

        if (loadViewMatrixAndZoom(viewMatrixFile, zoom, matrixString)) {
            // Set zoom and matrix text into UI
            w.findChild<QLineEdit*>("zoomValue")->setText(QString::number(zoom, 'f', 4));
            w.findChild<QTextEdit*>("matrixValue")->setPlainText(matrixString);

            // Trigger your existing update slots
            w.onZoomValueChanged();
            w.onMatrixValueChanged();
        }
    }

    w.show();

    if (qApp->arguments().indexOf(QRegExp("-writefib")) != -1 ||
        qApp->arguments().indexOf(QRegExp("-writeobj")) != -1 ||
        qApp->arguments().indexOf(QRegExp("-screenshot")) != -1) {
        return 0;
    }

    return a.exec();
}
