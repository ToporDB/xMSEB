#include <QtCore/QCoreApplication>
#include "connections.h"
//#include <qfile.h>
//#include <qtextstream.h>
#include <QException>
#include <QtDebug>

#include <QStringList>

#include "qmath.h"

#include <QTextStream>
#include <QFile>
#include <QDataStream>

Connections::~Connections() {
    // Delete all dynamically allocated Edge objects
    for (Edge* edge : qAsConst(edges)) {
        delete edge;
    }

    // Free allocated memory for comps and directions
    delete[] comps;
    delete[] directions;
}

Connections::Connections(QString nname, QString ename, QString fileName)
{
    params();
    prefix = fileName; //nname;
    QFile n(nname);
    qDebug() << nname;
    if (!n.open(QIODevice::ReadOnly)) qDebug("nodes unreadable");
    QTextStream ns(&n);
    QString nl;

    while(!ns.atEnd()) {
        nl = ns.readLine();

        QStringList vals = nl.split(" ", Qt::SkipEmptyParts);
        QVector3D* anode;
        //x,y,z
        anode = new QVector3D(((QString)(vals.at(0))).toDouble(),
                              ((QString)(vals.at(1))).toDouble(),
                              ((QString)(vals.at(2))).toDouble());
        // qDebug() << anode->x() << anode->y() << anode->z();
        nodes << *anode;

    }
    n.close();
    qDebug() << "nodes read";

    QFile e(ename);
    if (!e.open(QIODevice::ReadOnly)) {
        qDebug("edges unreadable");
        exit(2);
    }
    QTextStream es(&e);
    QString el;
    while (!es.atEnd()) {
        el = es.readLine();
        QStringList evals = el.split(" ", Qt::SkipEmptyParts);

        if (evals.size() < 2) {
            qDebug() << "Invalid edge line: " << el;
            continue;
        }

        int from = evals.at(0).toInt();
        int to = evals.at(1).toInt();

        const double defaultWeight = 1.0;
        double weight = (evals.size() > 2) ? evals.at(2).toDouble() : defaultWeight;

        // Check if start and end clusters are given; otherwise, set to "dummy"
        QString startCluster = (evals.size() > 3) ? evals.at(3) : "dummy";
        QString endCluster = (evals.size() > 4) ? evals.at(4) : "dummy";

        // Negative weights
        if (weight < 0.0) {
            std::swap(from, to);
            std::swap(startCluster, endCluster);
            weight = std::abs(weight);
        }
        // Adding 1 so we don't divide with number between (0, 1)
        weight += 1.0;

        try {
            Edge* aedge = new Edge(nodes.at(from), nodes.at(to), QString::number(weight), startCluster, endCluster);
            edges << aedge;
        } catch (QException& e) {
            qDebug() << "Out of bounds error for nodes:" << from << to;
        }
    }
    e.close();

    qDebug() << edges.length() << " edges...";
}

Connections::Connections(QString fib){
    //TODO: Das hier so umbiegen, dass ich Gabys Daten laden kann...
    params();
    prefix = fib;
    QFile n(fib);

    if (!n.open(QIODevice::ReadOnly)) {
        qDebug() << "vtk unreadable: " << fib;
        exit(1);
    }

    QTextStream ns(&n);
    QString nl;
    QDataStream ins(&n);
    ins.setByteOrder(QDataStream::BigEndian);
    ins.setFloatingPointPrecision(QDataStream::SinglePrecision);
    nl = ns.readLine(); //skip first lines;
    nl = ns.readLine(); //TODO: Other types of stuff...
    nl = ns.readLine();
    nl = ns.readLine();
    nl = ns.readLine();

    //ns.pos();
    qDebug() << ns.pos(); //TODO: Das hier sollte nichts ausmachen, tuts aber...
    qDebug() << nl;
    QStringList vals = nl.split(" ");
    float np = ((QString)(vals.at(1))).toInt();
    qDebug() << np;
    for (int i = 0; i < np; i++) {
        QVector3D* anode;
        float x,y,z;
        ins >> x;
        ins >> y;
        ins >> z;
        anode = new QVector3D(x,y,z);
        nodes << *anode;
    }
    qDebug() << ns.pos(); //TODO: WTF, siehe oben?
    ns.seek(n.pos() + np*3*4 + 1); //Textstream aufs Zeichen nach den Punkten...
    qDebug() << ns.pos();
    nl = ns.readLine();
    qDebug() << nl;
    vals = nl.split(" ");
    float ncons = ((QString)(vals.at(1))).toInt();

    qDebug() << ns.pos();
    for (int i = 0; i < ncons; i++) {
        qint32 numpoints;
        QString weight;
        weight = "10";

        QString startCluster;
        startCluster = "dummy";

        QString endCluster;
        endCluster = "dummy";
        ins >> numpoints;
        //qDebug() << numpoints;
        qint32* ps = new qint32[numpoints];
        for (int pn = 0; pn < numpoints; pn++){
            ins >> ps[pn];
        }
        Edge* aedge = new Edge(nodes.at(ps[0]), nodes.at(ps[numpoints-1]), weight, startCluster, endCluster);
        aedge->points.removeLast();
        for (int pn = 1; pn < numpoints; pn++){
            aedge->points << nodes.at(ps[pn]);
        }
        edges << aedge;
    }
    n.close();
    qDebug() << "nodes read";
}

void Connections::params() {
    c_thr = 0.8;
    start_i = 10;
    numcycles = 10;
    bell = 5;
    smooth = 3;
    beta = 0.1;
    checkpoints = 0;
    lane_width = 1.0f;
    directed = 0;
}

void Connections::subdivide(int newp) {
    for (int i = 0; i < edges.size(); ++i) {
        Edge* e = edges.at(i);

        if (e->points.size() < newp) {
            e->subdivide(newp);
        }
    }
}

double Connections::attract() {
    double totalMovement = 0.0;

    for (int ie = 0; ie < edges.size(); ++ie) {
        Edge* e = edges.at(ie);
        double weightOfThisEdge = e->wt.toDouble();
        QVector<QVector3D> newForces(e->points.length());  // Thread-local forces

        double edgeMovement = 0.0;

#pragma omp parallel for reduction(+:edgeMovement)
        for (int i = 1; i < e->points.length() - 1; ++i) {
            QVector3D p = e->points.at(i);
            double edgeDepthFactor = (-i * (i - e->points.length())) / std::pow(e->points.length() / 2.0, 2);
            double fsum = 0;
            QVector3D f(0, 0, 0);

            QVector3D e_dir = e->points.last() - e->points.first();
            e_dir.normalize();

            for (int ef = 0; ef < edges.size(); ++ef) {
                float c = comp(ie, ef);
                if (c <= c_thr || (directed && directions[ie + edges.size() * ef] < 0)) continue;

                Edge* other = edges.at(ef);

                QVector3D potential;
                double weight = 0.0;

                std::tie(potential, weight) = computeUndirectedAttractionForce(e, other, i, ie, ef);

                fsum += weight;
                f += weight * potential;
            }

            f /= fsum;
            QVector3D force = edgeDepthFactor * edgeDepthFactor * ((f - p) / weightOfThisEdge);
            force += beta * e->forces.at(i);  // Momentum
            newForces[i] = force;

            edgeMovement += force.length();
        }

        // Apply the computed forces sequentially
        for (int i = 1; i < e->points.length() - 1; ++i) {
            e->forces[i] = newForces[i];
        }

        totalMovement += edgeMovement;
    }

    // Safe to apply forces now
    for (Edge* e : qAsConst(edges)) {
        e->applyForces();
    }

    return totalMovement;
}


void Connections::fullAttract() {
    calcComps();

    double spfac = 1.3;
    double spnow = 1;
    int i = start_i;
    for (int cycle = 0; cycle < numcycles; cycle++){
        int sps = qRound(spnow);
        subdivide(sps);
        qDebug() << "starting " << i << " iterations with c_thr:" << c_thr << "segments: " << edges.first()->points.length()-1;
        for (int j = 0; j<i; j++){
            double movement = attract();
            if (checkpoints) writeVTK(j, cycle);

            if (movement < lambda && qRound(spnow) != 1) {
                qDebug() << "Layout converged after" << cycle << "cycles and" << j << "iterations.";
                break;
            }
        }
        i = std::max(1, i - 1);
        spnow *= spfac;
    }

    if (directed) {
        addLateralForce();

        for (int i = 0; i < start_i; ++i) {
            attract();
        }
    }

    // for further subdivision without attraction
    if (numcycles > 0 && smooth > 1) {
        for (int i = 0; i < smooth; ++i) {
            subdivide(qRound(spnow) + i);
            qDebug() << "Number of subdivision points:" << qRound(spnow) + i;
        }
    }
}

void Connections::addLateralForce() {
    for (int ei = 0; ei < edges.length(); ++ei) {
        Edge* p = edges.at(ei);
        double weightOfThisEdge = p->wt.toDouble();
        QVector<QVector3D> newForces(p->points.length());


#pragma omp parallel for
        for (int i = 1; i < p->points.length() - 1; ++i) {
            QVector3D p_i = p->points.at(i);
            double edgeDepthFactor = (-i * (i - p->points.length())) / std::pow(p->points.length() / 2.0, 2);
            double fsum = 0;
            QVector3D f(0, 0, 0);

            for (int ej = ei + 1; ej < edges.length(); ++ej) {
                float c = comp(ei, ej);
                if (c < c_thr) continue;

                Edge* q = edges.at(ej);

                QVector3D potential;
                double weight = 0.0;

                std::tie(potential, weight) = computeDirectedAttractionForce(p, q, i, ei, ej);

                fsum += weight;
                f += weight * potential;
            }

            if (fsum > 0) {
                f /= fsum;
                QVector3D force = edgeDepthFactor * edgeDepthFactor * ((f - p_i) / weightOfThisEdge);
                newForces[i] = force;
            }
        }

        // Apply the computed forces sequentially
        for (int i = 1; i < p->points.length() - 1; ++i) {
            p->forces[i] = newForces[i];
        }
    }

    for (Edge* e : qAsConst(edges)) {
        e->applyForces();
    }
}

void Connections::calcComps(){
    comps = new float[edges.size()*edges.size()];
    if (directed) directions = new float[edges.size()*edges.size()];

    #pragma omp parallel for num_threads (7)
    for (int i=0; i<edges.length(); i++){
        for (int j=0; j<edges.length(); j++){
            if (i==j) {
                comps[i+edges.size()*j]=0.0;
            } else {
                Edge* ei = edges.at(i);
                Edge* ej = edges.at(j);

                //calculate compatibility btw. edge i and j
                //angle
                double angle_comp;
                if (!ei->flip(ej)) {
                    angle_comp = QVector3D::dotProduct(ei->fn-ei->tn,ej->tn-ej->fn);
                } else {
                    angle_comp = QVector3D::dotProduct(ei->fn-ei->tn,ej->fn-ej->tn);
                }
                angle_comp /= ei->length()*ej->length();

                // New: Precalculate the directionality of two edges, so we can use it for diretionality if needed
                if (directed) directions[i + edges.size() * j] = QVector3D::dotProduct(
                        (ei->points.last() - ei->points.first()).normalized(),
                        (ej->points.last() - ej->points.first()).normalized()
                    );

                //length
                double lavg = (ei->length()+ej->length())/2.0;
                double l_comp = 2 / ((lavg/qMin(ei->length(),ej->length())) + (qMax(ei->length(),ej->length())/lavg));
                //position
                QVector3D mi = (ei->fn+ei->tn)/2;
                QVector3D mj = (ej->fn+ej->tn)/2;
                double p_comp = lavg / (lavg + (mi-mj).length());
                //visibility
                if (angle_comp * l_comp * p_comp > 0.9) {
                    double vis_comp = qMin(vis_c(ei,ej),vis_c(ej,ei));
                    comps[i+edges.size()*j] = angle_comp * l_comp * p_comp * vis_comp;
                } else {
                    comps[i+edges.size()*j] = angle_comp * l_comp * p_comp;
                }

            }
        }
    }
}

double Connections::vis_c(Edge* ep, Edge* eq) {
    QVector3D i0 = proj(ep->fn,ep->tn, eq->fn);
    QVector3D i1 = proj(ep->fn,ep->tn, eq->tn);
    QVector3D im = (i0+i1)/2;
    QVector3D pm = (ep->fn+ep->tn)/2;

    return qMax(1-2*(pm-im).length()/(i0-i1).length(),(float)0.0);
}

QVector3D Connections::proj(QVector3D a, QVector3D b, QVector3D p) {
    QVector3D ba = b-a;
    QVector3D pa = p-a;
    return a + ba*QVector3D::dotProduct(ba,pa) / ba.lengthSquared();
}

float Connections::comp(int i, int j) const {
    return comps[i+edges.size()*j];
}

void Connections::writeVTK() {
    writeVTK(-1, -1);
}

void Connections::writeVTK(int current_start_i = -1, int current_numcycle = -1) {
    qDebug() << "Writing VTK file...";

    QFile file(name(current_start_i, current_numcycle) + ".vtk");
    if (!file.open(QIODevice::WriteOnly)) {
        qWarning() << "Error opening file for writing:" << file.fileName();
        return;
    }
    QTextStream out(&file);
    out.setRealNumberPrecision(10);

    int totalPoints = 0;
    int totalLines = edges.size();

    // First, calculate the total number of points dynamically
    for (const Edge* edge : qAsConst(edges)) {
        totalPoints += edge->points.size();
    }

    // Write VTK Header
    out << "# vtk DataFile Version 3.0" << Qt::endl;
    // out << "Generated by Connections::writeVTK()" << Qt::endl;
    out << "ASCII" << Qt::endl;
    out << "DATASET POLYDATA" << Qt::endl;
    out << "POINTS " << totalPoints << " float" << Qt::endl;

    // Write all points
    int pointIndex = 0;  // Track global point index
    for (const Edge* edge : qAsConst(edges)) {
        for (const QVector3D& po : edge->points) {
            out << po.x() << " " << po.y() << " " << po.z() << Qt::endl;
            ++pointIndex;
        }
    }

    // Write the connectivity (LINES)
    out << "LINES " << totalLines << " " << (totalPoints + totalLines) << Qt::endl;
    pointIndex = 0;
    for (const Edge* edge : qAsConst(edges)) {
        int edgeSize = edge->points.size();
        if (edgeSize == 0) continue;  // Skip empty edges

        out << edgeSize;
        for (int i = 0; i < edgeSize; i++) {
            out << " " << pointIndex++;
        }
        out << Qt::endl;
    }

    file.close();
    qDebug() << "VTK file successfully written:" << file.fileName();
}

void Connections::writeBinaryVTK(){
    qDebug() << "writing binary vtk file";
    writeBinaryVTK(name());
}

void Connections::writeBinaryVTK(QString name){

    qDebug() << "writing: " << name;
    QFile file(name+".vtk");
    if (!file.open(QIODevice::WriteOnly)) qDebug() << "error opening file for writing";
    QDataStream out(&file);
    QTextStream outt(&file);

    out.setByteOrder(QDataStream::BigEndian);
    out.setFloatingPointPrecision(QDataStream::SinglePrecision);

    int n = edges.size();
    int m = edges.at(0)->points.size();

    outt << "# vtk DataFile Version 3.0" << Qt::endl;
    outt << "I am a header! Yay!" << Qt::endl;
    outt << "BINARY" << Qt::endl;
    outt << "DATASET POLYDATA" << Qt::endl;
    outt << "POINTS " << m*n << " float" << Qt::endl;

    for (int e = 0; e<n; e++){
        Edge* ed = edges.at(e);
        for (int p=0; p<ed->points.size(); p++){
            QVector3D po = ed->points.at(p);
            out << (float)po.x() << (float)po.y() << (float)po.z();
        }
    }
    outt << Qt::endl;

    outt << "LINES " << n << " " << n*(m+1) << Qt::endl;
    int i = 0;
    for (int e = 0; e<n; e++){
        out << m;
        for (int p=0; p<m; p++){
            out << i++;
        }
    }
    outt << Qt::endl;

    file.close();

}


void Connections::writeSegments(){
    qDebug() << "writing segments file";

    QFile file("segments");
    file.open(QIODevice::WriteOnly);
    QTextStream out(&file);

    int n = edges.size();

    for (int e = 0; e<n; e++){
        Edge* ed = edges.at(e);
        for (int p=0; p<ed->points.size(); p++){
            QVector3D po = ed->points.at(p);
            out << (float)po.x() << " " << (float)po.y()  << " " << (float)po.z() << " " << e << Qt::endl;
        }
    }

    file.close();
    qDebug() << "file written";
}

QString Connections::name() {
    return name(-1, -1);
}

QString Connections::name(int current_start_i = -1, int current_numcycles = -1) {
    int output_start_i = current_start_i != -1 ? current_start_i : start_i;
    int output_numcycles = current_numcycles != -1 ? current_numcycles : numcycles;
    return prefix +
           "_c_thr" + QString::number(c_thr,'f',4) +
           "_numcycles" + QString("%1").arg(output_numcycles,2,10,QLatin1Char('0') ) +
           "_start_i" + QString("%1").arg(output_start_i,4,10,QLatin1Char('0'))+
           "_directed" + QString::number(directed);
}

std::pair<QVector3D, double> Connections::computeUndirectedAttractionForce(
    Edge* e, Edge* other, int i, int ei, int ej
    ) const {
    QVector3D p = e->points.at(i);
    QVector3D pe = other->flip(e)
                       ? other->points.at(i)
                       : other->points.at(other->points.length() - 1 - i);

    double weight = qExp(-(pe - p).lengthSquared() / (2 * bell * bell)) /
                    other->wt.toDouble() * std::pow(comp(ei, ej), 2);

    return {pe, weight};
}

std::pair<QVector3D, double> Connections::computeDirectedAttractionForce(
    Edge* e, Edge* other, int i, int ei, int ej
    ) const {
    float directionFactor = directions[ei + edges.length() * ej];

    // Choose index based on direction
    int other_i = (directionFactor < 0)
                      ? other->points.length() - 1 - i
                      : i;

    QVector3D q_j = other->points.at(other_i);
    QVector3D p_i = e->points.at(i);

    if (directionFactor >= 0 || ((p_i - q_j).length() > lane_width * ((e->length() + other->length()) / 2) / 50)) {
        return {QVector3D(), 0.0f};
    }

    QVector3D q_prev = other->points.at(other_i - 1);
    QVector3D q_next = other->points.at(other_i + 1);

    QVector3D Tj = q_next - q_prev;
    if (Tj.lengthSquared() < 1e-6f)
        Tj = (other->points.last() - other->points.first());

    Tj.normalize();

    QVector3D up(0, 1, 0);
    QVector3D Nj = QVector3D::crossProduct(Tj, up).normalized();

    QVector3D potential = p_i + Nj * lane_width * (((e->length() + other->length()) / 2) / 50);

    double weight = qExp(-(potential - p_i).lengthSquared() / (2 * bell * bell)) /
                    other->wt.toDouble() * std::pow(comp(ei, ej), 2);

    return {potential, weight};
}
