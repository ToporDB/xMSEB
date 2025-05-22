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
    ~Connections();
    void params();
    void subdivide();
    void subdivide(int newp);
    double attract();
    void addLateralForce();
    void fullAttract();
    void calcComps();
    float* comps;
    float* directions;
    float comp(int i, int j) const;
    void writeVTK();
    void writeVTK(int, int);
    void writeBinaryVTK();
    void writeBinaryVTK(QString name);
    void writeSegments();
    QString name();
    QString name(int, int);

    // Computes attraction force for undirected edges
    std::pair<QVector3D, double> computeUndirectedAttractionForce(
        Edge* e, Edge* other, int i, int ei, int ej
        ) const;

    // Computes attraction force for directed edges
    std::pair<QVector3D, double> computeDirectedAttractionForce(
        Edge* e, Edge* other, int i, int ei, int ej
        ) const;

    double c_thr, bell, beta, lane_width, lambda = 1e-4;
    int start_i, numcycles, smooth, checkpoints, directed;
    QString prefix;

private:
    QList<QVector3D> nodes;
    QList<Edge*> edges;
    double vis_c(Edge* e1, Edge* e2);
    QVector3D proj(QVector3D a, QVector3D b, QVector3D p);
};

#endif // CONNECTIONS_H
