#include <QtCore/QCoreApplication>

#include <QtDebug>
#include <QStringList>
#include "connections.h"
#include <chrono>

QString arg(QString argname) {
    int nodespos = qApp->arguments().indexOf(QRegExp("-"+argname+"*"));
    if (nodespos!=-1) {
        return qApp->arguments().at(nodespos+1);
    } else {
        return "";
    }
}

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();
    QCoreApplication a(argc, argv);

    if (qApp->arguments().indexOf(QRegExp("-help"))!=-1 || qApp->arguments().length()==1) {
            qDebug() << "bundler (-nodes <nodes> -cons <connections> -fileName <fileName> / -fib <.fib-file>) "
                    "[-c_thr <compatibility threshold>] "
                    "[-start_i <iterations in 1st cycle>] "
                    "[-numcycles <number of cycles>] "
                    "[-checkpoints <0 or 1>, default is 0] "
                    "[-directed <0 or 1>, default is 0] "
                    "[-bundles <0 or 1>, default is 0]";
            return 1;
    }

    qDebug() << "nodesname: " + arg("nodes");
    qDebug() << "consname: " + arg("cons");
    qDebug() << "filename: " + arg("fileName");

    Connections* cons;
    if (arg("cons")!="") {
        cons = new Connections(arg("nodes"),arg("cons"), arg("fileName"));
    } else if (arg("fib")!="") {
        cons = new Connections(arg("fib"));
    } else {
        qDebug() << "no dataset given";
        return 0;
    }

    if (arg("c_thr")!="") cons->c_thr = arg("c_thr").toDouble();
    if (arg("start_i")!="") cons->start_i = arg("start_i").toInt();
    if (arg("numcycles")!="") cons->numcycles = arg("numcycles").toInt();
    if (arg("bell")!="") cons->bell = arg("bell").toDouble();
    if (arg("smooth")!="") cons->smooth = arg("smooth").toInt();
    if (arg("checkpoints")!="") cons->checkpoints = arg("checkpoints").toInt();
    if (arg("directed")!="") cons->directed = arg("directed").toInt();
    if (arg("bundles")!="") cons->bundles = arg("bundles").toInt();
    if (arg("poly_deg")!="") cons->poly_deg = std::clamp(arg("poly_deg").toDouble(), 0.0, 3.0);

    qDebug() << cons->name();

    cons->fullAttract();
    cons->writeVTK();

    if (cons->bundles) cons->writeBundles();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    qDebug() << "Execution time: " << duration.count() << " ms";
    return 1;
}
