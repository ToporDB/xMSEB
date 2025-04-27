#ifndef CONNECTIONS_H
#define CONNECTIONS_H

#include <QString>
#include <QList>
#include "edge.h"
#include <QVector3D>

class Connections
{
public:
    Connections(QString nname, QString ename, QString fileName);
    Connections(QString fib);
    void params();
    void subdivide();
    void subdivide(int newp);
    void attract();

    void fullAttract();
    void calcComps();
    float* comps;
    float comp(int i, int j);
    void writeVTK();
    void writeVTK(int, int);
    void writeBinaryVTK();
    void writeBinaryVTK(QString name);
    void writeSegments();
    QString name();
    QString name(int, int);

    double c_thr, bell, beta;
    int start_i, numcycles, smooth;
    QString prefix;

private:
    QList<QVector3D> nodes;
    QList<Edge*> edges;
    double vis_c(Edge* e1, Edge* e2);
    QVector3D proj(QVector3D a, QVector3D b, QVector3D p);
};

#endif // CONNECTIONS_H
