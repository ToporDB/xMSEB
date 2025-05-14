#ifndef SPLINEUTILS_H
#define SPLINEUTILS_H

#include <QVector3D>
#include <QList>

namespace SplineUtils {

// Computes Catmull-Rom spline interpolated points
QList<QVector3D> interpolateCatmullRom(const QList<QVector3D>& controlPoints, int samplesPerSegment);

} // namespace SplineUtils
#endif // SPLINEUTILS_H
