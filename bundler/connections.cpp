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
    for (Edge* edge : qAsConst(edges)) {
        delete edge;
    }

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
    QStringList longName = fib.split('/');

    prefix = longName[longName.length() - 1];
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
    // TODO: set beta to 0.65f
    beta = 0.0;
    checkpoints = 0;
    lane_width = 1.0f;
    directed = 0;
    lambda = 1e-4;
    bundles = 0;
    poly_deg = 0.0;
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
    #pragma omp parallel for reduction(+:totalMovement)
    for (int ie = 0; ie < edges.size(); ++ie) {
        Edge* e = edges.at(ie);
        double weightOfThisEdge = e->wt.toDouble();

        for (int i = 1; i < e->points.length() - 1; ++i) {
            QVector3D p = e->points.at(i);
            double edgeDepthFactor = std::pow((4 * i * (e->points.length() - i)) / std::pow(e->points.length(), 2), poly_deg);
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

                std::tie(potential, weight) = computeUndirectedAttractionForce(e, other, i);

                fsum += weight;
                f += weight * potential;
            }

            if (fsum > 0) {
                f /= fsum;
                QVector3D force = edgeDepthFactor * edgeDepthFactor * ((f - p) / weightOfThisEdge);
                force += beta * e->forces.at(i);  // Momentum
                e->forces[i] = force;
                totalMovement += force.length();
            }
        }
    }

    for (Edge* e : qAsConst(edges)) {
        e->applyForces();
    }

    return totalMovement;
}

void Connections::fullAttract() {
    calcComps();

    double spfac = 1.3;
    double spnow = 2.028;
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
            direction_switch = !direction_switch;
        }
        i = std::max(1, i - 1);
        spnow *= spfac;
    }

    if (directed) {
        addLateralForces();
        attract();
    }

    // for further subdivision without attraction
    if (numcycles > 0 && smooth > 1) {
        for (int i = 0; i < smooth; ++i) {
            subdivide(qRound(spnow) + i);
            qDebug() << "Number of subdivision points:" << qRound(spnow) + i;
        }
    }
}

void Connections::addLateralForces() {
#pragma omp parallel for
    for (int ei = 0; ei < edges.length(); ++ei) {
        Edge* p = edges.at(ei);
        double weightOfThisEdge = p->wt.toDouble();


        for (int i = 1; i < p->points.length() - 1; ++i) {
            QVector3D p_i = p->points.at(i);
            double edgeDepthFactor = (-i * (i - p->points.length())) / std::pow(p->points.length() / 2.0, 2);
            double fsum = 0;
            QVector3D f(0, 0, 0);

            for (int ej = ei + 1; ej < edges.length(); ++ej) {
                float c = comp(ei, ej);

                if (c < c_thr || directions[ei + edges.size() * ej] >= 0) continue;

                Edge* q = edges.at(ej);

                QVector3D potential;
                double weight = 0.0;

                std::tie(potential, weight) = computeDirectedAttractionForce(p, q, i);

                fsum += weight;
                f += weight * potential;
            }

            if (fsum > 0) {
                f /= fsum;
                QVector3D force = edgeDepthFactor * edgeDepthFactor * ((f - p_i) / weightOfThisEdge);
                p->forces[i] = force;
            }
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

QVector3D Connections::proj(QVector3D a, QVector3D b, QVector3D p) const {
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

    for (const Edge* edge : qAsConst(edges)) {
        totalPoints += edge->points.size();
    }

    // Write VTK Header
    out << "# vtk DataFile Version 3.0" << Qt::endl;
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
        if (edgeSize == 0) continue;

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

void Connections::writeBundles() {
    qDebug() << "writing edge bundles file";

    QFile file(prefix + "_edge_bundles.txt");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qWarning() << "Failed to open output file";
        return;
    }

    QTextStream out(&file);

    int n = edges.size();
    float threshold = 0.1f;

    // Track which edges are already part of a bundle
    QVector<bool> edgeBundled(n, false);
    int bundle_id = 0;

    auto sign = [](float x) { return (x >= 0) ? 1 : -1; };

    for (int i = 0; i < n; ++i) {
        if (edgeBundled[i]) continue;

        QList<int> bundle_edges;
        float total_weight = 0.0f;
        float dir_sign = 1.0f;

        if (directed) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                float c = comp(i, j);
                if (c >= threshold) {
                    dir_sign = sign(directions[i * n + j]);
                    break;
                }
            }
        }

        for (int k = 0; k < n; ++k) {
            if (edgeBundled[k]) continue;
            float ck = comp(i, k);
            float dk = directed ? directions[i * n + k] : 1.0f;

            if (ck >= threshold && (!directed || sign(dk) == dir_sign)) {
                bundle_edges.append(k);
                edgeBundled[k] = true;
                total_weight += ck;
            }
        }

        if (bundle_edges.size() <= 1) continue;

        QString label = (dir_sign > 0) ? "forward" : "reverse";
        out << "Bundle " << bundle_id++ << " (" << label << "):\n";
        out << "  Edges: ";
        for (int idx : bundle_edges)
            out << idx << " ";
        out << "\n";
        out << "  Total Weight: " << total_weight << "\n\n";
    }

    file.close();
    qDebug() << "file written";
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
           "_numcycles" + QString("%1").arg(output_numcycles,2,10,QLatin1Char('0')) +
           "_start_i" + QString("%1").arg(output_start_i,4,10,QLatin1Char('0'))+
           "_directed" + QString::number(directed);
}

std::pair<QVector3D, double> Connections::computeUndirectedAttractionForce(
    Edge* e, Edge* other, int& i
    ) const {
    const QVector3D& p = e->points.at(i);

    int best_j = closest_intermediatePoint_index(other, e, i);
    if (best_j == -1) return {{0, 0, 0}, 0.0};  // no visible match

    const QVector3D& q = other->points.at(best_j);

    double weight = qExp(-(q - p).lengthSquared() / (2 * bell * bell)) / other->wt.toDouble();
    return {q, weight};
}


std::pair<QVector3D, double> Connections::computeDirectedAttractionForce(
    Edge* e, Edge* other, int &i
    ) const {
    // int other_i = other->points.length() - 1 - i;
    int other_i = closest_intermediatePoint_index(other, e, i);
    if (other_i == -1) return {{0, 0, 0}, 0.0};  // no visible match

    QVector3D q_j = other->points.at(other_i);
    QVector3D p_i = e->points.at(i);

    if (((p_i - q_j).length() > lane_width * ((e->length() + other->length()) / 2) / 50)) {
        return {QVector3D(), 0.0f};
    }

    QVector3D q_prev = other->points.at(other_i + 1);
    QVector3D q_next = other->points.at(other_i - 1);

    QVector3D T_j = q_next - q_prev;
    if (T_j.lengthSquared() < 1e-6f)
        T_j = (other->points.last() - other->points.first());

    T_j.normalize();

    QVector3D up(0, 1, 0);
    QVector3D N_j = QVector3D::crossProduct(T_j, up).normalized();

    QVector3D potential = q_j + N_j * lane_width * (((e->length() + other->length()) / 2) / 50);

    double weight = qExp(-(potential - p_i).lengthSquared() / (2 * bell * bell)) /
                    other->wt.toDouble();

    return {potential, weight};
}

bool Connections::vis_point_on_edge(Edge* Q, const QVector3D& p_i) const {
    const QVector3D& a = Q->points.first();
    const QVector3D& b = Q->points.last();

    QVector3D projPoint = proj(a, b, p_i);

    // Check if projection lies between a and b using dot products
    QVector3D ab = b - a;
    QVector3D ap = projPoint - a;
    QVector3D bp = projPoint - b;

    // If projection lies "between" a and b, both dot products should be ≤ 0
    return QVector3D::dotProduct(ab, bp) <= 0 && QVector3D::dotProduct(-ab, ap) <= 0;
}

int Connections::closest_intermediatePoint_index(Edge* Q, Edge* P, const int& i) const {
    const QVector3D& p = P->points.at(i);

    // Try matching index if it is "visible"
    bool flipOther = Q->flip(P);
    int j_match = flipOther ? i : Q->points.length() - 1 - i;

    if (j_match > 0 && j_match < Q->points.length() - 1) {
        const QVector3D& q_match = Q->points.at(j_match);
        bool p_visible_to_Q = vis_point_on_edge(Q, p);
        bool q_visible_to_P = vis_point_on_edge(P, q_match);

        if (p_visible_to_Q && q_visible_to_P) {
            return j_match;
        }
    }

    // If matching index is not visible, search closest visible point (excluding endpoints)
    int closest_j = -1;
    double minDistSq = std::numeric_limits<double>::max();

    // To have the best bundling result and not be biased by going only one way,
    // We are going to alter the direction each cycle.
    if (direction_switch) {
        for (int j = 1; j < Q->points.length() - 1; ++j) {
            const QVector3D& candidate = Q->points.at(j);
            if (!vis_point_on_edge(P, candidate)) continue;

            double d2 = (candidate - p).lengthSquared();
            if (d2 < minDistSq) {
                minDistSq = d2;
                closest_j = j;
            }
        }
    }
    else
    {
        for (int j = Q->points.length() - 1; j > 0 ; --j) {
            const QVector3D& candidate = Q->points.at(j);
            if (!vis_point_on_edge(P, candidate)) continue;

            double d2 = (candidate - p).lengthSquared();
            if (d2 < minDistSq) {
                minDistSq = d2;
                closest_j = j;
            }
        }
    }

    return closest_j;
}
