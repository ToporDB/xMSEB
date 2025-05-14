#include "SplineUtils.h"

namespace SplineUtils {

// Interpolates between 4 points using Catmull-Rom spline
static QVector3D catmullRom(const QVector3D& p0, const QVector3D& p1, const QVector3D& p2, const QVector3D& p3, float t) {
    float t2 = t * t;
    float t3 = t2 * t;

    return 0.5f * ((2.0f * p1) +
                   (-p0 + p2) * t +
                   (2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t2 +
                   (-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t3);
}

QList<QVector3D> interpolateCatmullRom(const QList<QVector3D>& controlPoints, int samplesPerSegment) {
    QList<QVector3D> result;

    int n = controlPoints.size();
    if (n < 2) return result;

    // Add the first straight segment
    result.append(controlPoints[0]);
    result.append(controlPoints[1]);

    if (n >= 4) {
        // Interpolate Catmull-Rom between interior segments
        for (int i = 1; i < n - 2; ++i) {
            for (int j = 0; j < samplesPerSegment; ++j) {
                float t = float(j) / samplesPerSegment;
                QVector3D point = catmullRom(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t);
                result.append(point);
            }
        }

        // Add the last straight segment
        result.append(controlPoints[n - 2]);
        result.append(controlPoints[n - 1]);
    }

    return result;
}

} // namespace SplineUtils
