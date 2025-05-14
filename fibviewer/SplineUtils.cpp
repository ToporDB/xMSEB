#include "SplineUtils.h"
#include <QtMath>

namespace SplineUtils {


// Helper function for centripetal Catmull-Rom parameterization
static float getT(float t, float alpha, const QVector3D& p0, const QVector3D& p1) {
    QVector3D d = p1 - p0;
    float a = QVector3D::dotProduct(d, d);
    return qPow(a, 0.5f * alpha) + t;
}

// Centripetal Catmull-Rom interpolation between 4 points
static QVector3D catmullRomCentripetal(const QVector3D& p0, const QVector3D& p1, const QVector3D& p2, const QVector3D& p3, float tNorm, float alpha = 0.5f) {
    float t0 = 0.0f;
    float t1 = getT(t0, alpha, p0, p1);
    float t2 = getT(t1, alpha, p1, p2);
    float t3 = getT(t2, alpha, p2, p3);

    float t = t1 + tNorm * (t2 - t1);

    QVector3D A1 = (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1;
    QVector3D A2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2;
    QVector3D A3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3;

    QVector3D B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2;
    QVector3D B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3;

    QVector3D C = (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2;

    return C;
}

QList<QVector3D> interpolateCatmullRom(const QList<QVector3D>& controlPoints, int samplesPerSegment, float alpha) {
    QList<QVector3D> result;

    int n = controlPoints.size();
    if (n < 2) return result;

    // Add first straight segment
    result.append(controlPoints[0]);
    result.append(controlPoints[1]);

    if (n >= 4) {
        for (int i = 1; i < n - 2; ++i) {
            for (int j = 0; j < samplesPerSegment; ++j) {
                float t = float(j) / samplesPerSegment;
                QVector3D point = catmullRomCentripetal(controlPoints[i - 1], controlPoints[i], controlPoints[i + 1], controlPoints[i + 2], t, alpha);
                result.append(point);
            }
        }
    }

    // Add last straight segment
    result.append(controlPoints[n - 2]);
    result.append(controlPoints[n - 1]);

    return result;
}


} // namespace SplineUtils
