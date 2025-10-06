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
    void addLateralForces();
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
    void writeBundles();
    QString name();
    QString name(int, int);
    bool vis_point_on_edge(Edge* Q, const QVector3D& p_i) const;
    int closest_intermediatePoint_index(Edge* Q, Edge* P, const int& i) const;
    std::pair<QVector3D, double> computeUndirectedAttractionForce(
        Edge* e, Edge* other, int &i
        ) const;

    std::pair<QVector3D, double> computeDirectedAttractionForce(
        Edge* e, Edge* other, int &i
        ) const;

    double c_thr, bell, beta, lane_width, lambda;
    int start_i, numcycles, smooth, checkpoints, directed, bundles, total_cycles;
    QString prefix;
    bool direction_switch = false;

private:
    QList<QVector3D> nodes;
    QList<Edge*> edges;
    double vis_c(Edge* e1, Edge* e2);
    QVector3D proj(QVector3D a, QVector3D b, QVector3D p) const;
};

#endif // CONNECTIONS_H
