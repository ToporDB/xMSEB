#include "connections.h"

#include <QtDebug>

#include <QStringList>

#include <QGLWidget>

Connections::Connections()
{
    min = QVector3D(-100,-100,-100);
    max = QVector3D(100,100,100);
    piv = (max-min)/2;
}

Connections::Connections(QString nname, QString ename)
{
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
        anode = new QVector3D(((QString)(vals.at(0))).toFloat(),
                              ((QString)(vals.at(1))).toFloat(),
                              ((QString)(vals.at(2))).toFloat());
        nodes << *anode;
    }
    n.close();
    qDebug() << "nodes read";

    QFile e(ename);
    if (!e.open(QIODevice::ReadOnly)) qDebug("edges unreadable");
    QTextStream es(&e);
    QString el;
    while(!es.atEnd()) {
        int f;
        int t;
        el = es.readLine();

        QStringList evals = el.split(" ", Qt::SkipEmptyParts);
        f = ((QString)(evals.at(0))).toInt();
        t = ((QString)(evals.at(1))).toInt();

        Edge* aedge;
        aedge = new Edge(nodes.at(f), nodes.at(t));
        if (aedge->length() > 0) edges << aedge;
    }
    e.close();

    qDebug() << edges.length() << " edges...";

    calculateBounds();
    createPrims();
}

Connections::Connections(QString fib) {
    QFile n(fib);
    if (!n.open(QIODevice::ReadOnly)) {
        qDebug() << "vtk unreadable: " << fib;
        return;
    }

    // Check if the file is ASCII or Binary by reading the header
    QTextStream ns(&n);
    QString nl;

    while (!ns.atEnd()) {
        nl = ns.readLine();
        if (nl.startsWith("ASCII")) {
            qDebug() << "Detected ASCII format";
            n.seek(0);
            readAscii(n);
            break;
        }
        if (nl.startsWith("BINARY")) {
            qDebug() << "Detected Binary format";
            n.seek(0);
            readBinary(n);
            break;
        }
    }

    n.close();
    calculateBounds();
    createPrims();
}

void Connections::readAscii(QFile& n) {
    QTextStream ns(&n);
    QString nl;

    while (!ns.atEnd()) {
        nl = ns.readLine();
        if (nl.startsWith("POINTS")) {
            break;
        }
    }

    QStringList vals = nl.split(" ");
    int np = vals[1].toInt();
    qDebug() << "Number of points: " << np;

    // Read points
    for (int i = 0; i < np; i++) {
        nl = ns.readLine();
        QStringList pointData = nl.split(" ", Qt::SkipEmptyParts);
        nodes.append(QVector3D(pointData[0].toFloat(), pointData[1].toFloat(), pointData[2].toFloat()));
    }

    // Read LINES section
    while (!ns.atEnd()) {
        nl = ns.readLine();
        if (nl.startsWith("LINES")) {
            break;
        }
    }

    vals = nl.split(" ");
    int ncons = vals[1].toInt();
    qDebug() << "Number of connections: " << ncons;

    // Read connections (LINES)
    for (int i = 0; i < ncons; i++) {
        nl = ns.readLine();
        QStringList connData = nl.split(" ", Qt::SkipEmptyParts);
        int numpoints = connData[0].toInt();

        Edge* aedge = new Edge(nodes.at(connData[1].toInt()), nodes.at(connData[numpoints].toInt()));
        aedge->points.removeLast();

        for (int j = 2; j <= numpoints; j++) {
            aedge->points << nodes.at(connData[j].toInt());
        }
        edges.append(aedge);
    }
}

void Connections::readBinary(QFile& n) {
    QString nl;
    QTextStream ns(&n);
    QDataStream ins(&n);
    ins.setByteOrder(QDataStream::BigEndian);
    ins.setFloatingPointPrecision(QDataStream::SinglePrecision);

    while (!ns.atEnd()) {
        nl = ns.readLine();
        if (nl.startsWith("POINTS")) {
            break;
        }
    }

    ns.pos();
    QStringList vals = nl.split(" ");
    int np = ((QString)(vals.at(1))).toInt();
    qDebug() << "number of points: " << np;
    for (int i = 0; i < np; i++) {
        QVector3D* anode;
        float x,y,z;
        ins >> x;
        ins >> y;
        ins >> z;
        anode = new QVector3D(x,y,z);
        nodes << *anode;
    }
    ns.pos();
    ns.seek(n.pos() + np*3*4 + 1);
    ns.pos();
    nl = ns.readLine();
    vals = nl.split(" ");
    int ncons = ((QString)(vals.at(1))).toInt();
    // int nps = ((QString)(vals.at(2))).toInt();
    qDebug() << "number of connections: " << ncons;

    ns.pos();
    for (int i = 0; i < ncons; i++) {
        qint32 numpoints;
        ins >> numpoints;
        qint32* ps = new qint32[numpoints];
        for (int pn = 0; pn < numpoints; pn++){
            ins >> ps[pn];
        }
        Edge* aedge = new Edge(nodes.at(ps[0]), nodes.at(ps[numpoints-1]));
        aedge->points.removeLast();
        for (int pn = 1; pn < numpoints; pn++){
            aedge->points << nodes.at(ps[pn]);
        }
        edges << aedge;
    }
}



void Connections::createPrims(){
    for (int i = 0; i < edges.size(); ++i) {
        Edge* e = edges.at(i);
        for (int i=0; i<e->points.length()-1; i++){
            QVector3D* p1 = new QVector3D(e->points.at(i));
            QVector3D* p2 = new QVector3D(e->points.at(i+1));

            Segment* seg = new Segment(p1,p2);
            prims << seg;
        }
    }
}

void Connections::calculateBounds() {
    if (edges.size() == 1) {
        piv = QVector3D(edges.at(0)->fn.x() + 2.0f, edges.at(0)->fn.y() + 2.0f, edges.at(0)->fn.z() + 2.0f);
        return;
    }

    for (int i = 0; i < edges.size(); i++) {
        Edge *aedge = edges.at(i);
        QVector3D afn = aedge->fn;
        QVector3D atn = aedge->tn;
        if (i==0){
            max = QVector3D(afn.x(),afn.y(),afn.z());
            min = QVector3D(afn.x(),afn.y(),afn.z());
        } else {
            max.setX(qMax(afn.x(),max.x()));
            max.setY(qMax(afn.y(),max.y()));
            max.setZ(qMax(afn.z(),max.z()));
            min.setX(qMin(afn.x(),min.x()));
            min.setY(qMin(afn.y(),min.y()));
            min.setZ(qMin(afn.z(),min.z()));

            max.setX(qMax(atn.x(),max.x()));
            max.setY(qMax(atn.y(),max.y()));
            max.setZ(qMax(atn.z(),max.z()));
            min.setX(qMin(atn.x(),min.x()));
            min.setY(qMin(atn.y(),min.y()));
            min.setZ(qMin(atn.z(),min.z()));
        }
    }
    piv = -(min+max)/2;
}

bool primLTprim(const Primitive* p1, const Primitive* p2) {
    return *p1 < *p2;
}

void Connections::paintPoints(){
    // if (nodes.isEmpty()) return;

    // // Draw first node (index 0): black, size 10
    // glColor3f(0.0f, 0.0f, 0.0f); // Black
    // glPointSize(10.0f);
    // glBegin(GL_POINTS);
    // QVector3D node = nodes.first();
    // glVertex3f(node.x(), node.y(), node.z());
    // node = nodes.last();
    // glVertex3f(node.x(), node.y(), node.z());
    // glEnd();

    // // Draw middle nodes: red, size 5
    // if (nodes.size() > 2) {
    //     glColor3f(1.0f, 0.0f, 0.0f); // Red
    //     glPointSize(5.0f);
    //     glBegin(GL_POINTS);
    //     for (int i = 1; i < nodes.size() - 1; ++i) {
    //         QVector3D n = nodes.at(i);
    //         glVertex3f(n.x(), n.y(), n.z());
    //     }
    //     glEnd();
    // }
}

void Connections::sortCons(){
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX,mat);
    for (int i = 0; i < prims.length(); i++){
        Segment* p = (Segment*)prims.at(i);
        double z1 = p->p1->x()*mat[2]+p->p1->y()*mat[6]+p->p1->z()*mat[10]+mat[14];
        double z2 = p->p2->x()*mat[2]+p->p2->y()*mat[6]+p->p2->z()*mat[10]+mat[14];
        p->depth = (z1+z2)/2.0;
    }
    std::sort(prims.begin(),prims.end(),primLTprim);
}

void Connections::writeBinaryVTK(QString filename){
    qDebug() << "writing binary vtk file";

    QFile file(filename);
    file.open(QIODevice::WriteOnly);
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
    qDebug() << "file written";
}

void Connections::writeCSVs(){
    QFile cfile("coords.csv");
    cfile.open(QIODevice::WriteOnly);
    QTextStream out(&cfile);
    out << "id,label,x,y" << Qt::endl;
    for (int i = 0; i<nodes.length(); i++){
        out << i << "," << i << "," << nodes.at(i).x() << "," << nodes.at(i).y() << Qt::endl;
    }
    cfile.close();

    QFile ffile("cons.csv");
    ffile.open(QIODevice::WriteOnly);
    QTextStream fout(&ffile);
    fout << "from,to,val" << Qt::endl;
    for (int i = 0; i<edges.length(); i++){
        fout << nodes.indexOf(edges.at(i)->fn) << "," << nodes.indexOf(edges.at(i)->tn) << "," << 1 << Qt::endl;
    }
    ffile.close();
}

void Connections::writeOBJ(QString filename){
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream out(&file);

    qDebug() << "writing vertices for " << edges.length() << " edges...";
    for (int e = 0; e<edges.length(); e++){
        Edge* ed = edges.at(e);
        qDebug() << e;
        for (int p=0; p<ed->points.size(); p++){
            QVector3D po = ed->points.at(p);
            out << "v " << (float)po.x() << " " <<(float)po.y() << " " << (float)po.z() << Qt::endl;
        }
    }

    //1st point is point 1 in .obj...
    int i = 1;
    qDebug() << "writing " << edges.length() << " edges...";
    for (int e = 0; e<edges.length(); e++){
        Edge* ed = edges.at(e);
        qDebug() << e;
        out << "l";
        for (int p=0; p<ed->points.size(); p++){
            out << " " << i;
            i++;
        }
        for (int p=0; p<ed->points.size()-2; p++){
            out << " " << (i-2)-p;
        }
        out << Qt::endl;
    }
    file.close();
}
