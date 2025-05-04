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
    QVector3D computeDirectionalPotential(
        const QVector3D& q_j,
        const QVector3D& q_prev,
        const QVector3D& q_next,
        const QVector3D& e_dir,
        const QVector3D& q_dir,
        float lane_width,
        double weightOfComparedEdge
        ) const;
    std::pair<QVector3D, double> computeDirectedAttractionForce(
        Edge* e, Edge* other, int i, const QVector3D& p,
        const QVector3D& e_dir, float c
        ) const;

    std::pair<QVector3D, double> computeUndirectedAttractionForce(
        Edge* e, Edge* other, int i, const QVector3D& p, float c
        ) const;

    double c_thr, bell, beta, lane_width;
    int start_i, numcycles, smooth, checkpoints, directed;
    QString prefix;

private:
    QList<QVector3D> nodes;
    QList<Edge*> edges;
    double vis_c(Edge* e1, Edge* e2);
    QVector3D proj(QVector3D a, QVector3D b, QVector3D p);
};

#endif // CONNECTIONS_H
