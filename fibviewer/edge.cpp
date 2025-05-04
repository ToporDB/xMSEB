#include "edge.h"

#include <QtDebug>

Edge::Edge(QVector3D fn, QVector3D tn)
{
    this->fn = fn;
    this->tn = tn;
    points << fn << tn;
}

void Edge::paintGL(bool intermediateNodes, bool startAndEndNodes, float alpha, bool dualGradient) {
    if (startAndEndNodes) {
        glPointSize(10);
        glColor4f(0.0, 0.0, 0.0, alpha);
        glBegin(GL_POINTS);
        glVertex(fn);
        glVertex(tn);
        glEnd();
    }

    glLineWidth(1);
    glPointSize(5);

    int n = points.length();
    if (n < 2) return;

    for (int i = 0; i < n - 1; ++i) {
        QVector3D p1 = points.at(i);
        QVector3D p2 = points.at(i + 1);

        QVector3D nor = p1 - p2;

        QVector3D col1, col2;

        if (dualGradient) {
            float t1 = float(i) / float(n - 1);
            float t2 = float(i + 1) / float(n - 1);

            // Red → Black → Blue gradient
            auto gradient = [](float t) -> QVector3D {
                if (t < 0.5f) {
                    float f = t / 0.5f; // 0 to 1
                    return QVector3D(1.0f * (1 - f), 0.0f, 0.0f); // Red to black
                } else {
                    float f = (t - 0.5f) / 0.5f; // 0 to 1
                    return QVector3D(0.0f, 0.0f, 1.0f * f); // Black to blue
                }
            };

            col1 = gradient(t1);
            col2 = gradient(t2);
        } else {
            // Solid red
            col1 = col2 = QVector3D(1.0f, 0.0f, 0.0f);
        }

        glBegin(GL_LINES);
        glNormal3f(nor.x(), nor.y(), nor.z());

        glColor4f(col1.x(), col1.y(), col1.z(), alpha);
        glVertex(p1);

        glColor4f(col2.x(), col2.y(), col2.z(), alpha);
        glVertex(p2);
        glEnd();

        if (intermediateNodes && i != 0) {
            glColor4f(col1.x(), col1.y(), col1.z(), alpha);
            glBegin(GL_POINTS);
            glVertex(p1);
            glEnd();
        }
    }
}

void Edge::glVertex(QVector3D v){
    glVertex3f(v.x(),v.y(),v.z());
}

double Edge::length(){
    return (fn-tn).length();
}

double Edge::minDist(QVector3D p1, QVector3D p2, QVector3D p){
    double l2 = (p1-p2).lengthSquared();
    if (l2 == 0.0) return (p-p1).length();
    QVector3D bla1 = p-p1;
    QVector3D bla2 = p2-p1;
    double t = QVector3D::dotProduct(bla1,bla2)/l2;
    if (t < 0.0) return (p1-p).length();
    else if (t > 1.0) return (p2-p).length();
    return (p-(p1+t*(p2-p1))).length();
}

double Edge::minDist(QVector3D p){
    double minD = 0.0;
    for (int i = 0; i<points.length()-1;i++){
        double distSeg = minDist(points.at(i),points.at(i+1),p);
        if (i==0) minD = distSeg;
        if (distSeg<minD) minD = distSeg;
    }
    return minD;
}
