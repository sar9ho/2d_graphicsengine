
/*
 *  Copyright 2025 Sarah Ouda
 */


#include "include/GMatrix.h"
#include <cmath>    // for sin, cos
#include <optional>
#include <cstring>

// initializeing identity matrix  [1, 0, 0]
                               // [0, 1, 0]
                               // [0, 0, 1]
GMatrix::GMatrix() {
    fMat[0] = 1; fMat[2] = 0; fMat[4] = 0;
    fMat[1] = 0; fMat[3] = 1; fMat[5] = 0;
}


//  translation matrix --- [1, 0, tx, 0, 1, ty]
// [1, 0, tx]
// [0, 1, ty]  tx and ty = how much to shift objects in x and y directions
// [0, 0, 1]
GMatrix GMatrix::Translate(float tx, float ty) {
    return GMatrix(1, 0, tx, 0, 1, ty);
}

// scaling matrix --- [sx, 0, 0, 0, sy, 0]
// [sx, 0, 0]
// [0, sy, 0]  sx scales along the x ax, and sy scales along the y ax
// [0, 0, 1]
GMatrix GMatrix::Scale(float sx, float sy) {
    return GMatrix(sx, 0, 0, 0, sy, 0);
}

// rotation matrix --- rotation by radians
// [cos(θ), -sin(θ), 0]
// [sin(θ), cos(θ), 0]  
// [0,       0,    1]
GMatrix GMatrix::Rotate(float radians) {
    float c = cos(radians);
    float s = sin(radians);
    return GMatrix(c, -s, 0, s, c, 0);
}

// multiply two matrices (a * b) /matrix mult
GMatrix GMatrix::Concat(const GMatrix& a, const GMatrix& b) {
    return GMatrix(
        a.fMat[0] * b.fMat[0] + a.fMat[2] * b.fMat[1],                    // new a
        a.fMat[0] * b.fMat[2] + a.fMat[2] * b.fMat[3],                    // new c
        a.fMat[0] * b.fMat[4] + a.fMat[2] * b.fMat[5] + a.fMat[4],          // new e
        a.fMat[1] * b.fMat[0] + a.fMat[3] * b.fMat[1],                    // new b
        a.fMat[1] * b.fMat[2] + a.fMat[3] * b.fMat[3],                    // new d
        a.fMat[1] * b.fMat[4] + a.fMat[3] * b.fMat[5] + a.fMat[5]           // new f
    );
}

// invert matrix -- computes inverse of the matrix using determinant
std::optional<GMatrix> GMatrix::invert() const {  //optional cus empty if not invert
    float det = fMat[0] * fMat[3] - fMat[1] * fMat[2];  // Compute determinant (a*d - b*c)
    if (det == 0) return {}; // No inverse exists.
    float invDet = 1.0f / det;
    return GMatrix(
        fMat[3] * invDet,                   // new a = d/det
        -fMat[2] * invDet,                  // new c = -c/det
        (fMat[2] * fMat[5] - fMat[3] * fMat[4]) * invDet,  // new e = (c*f - d*e)/det
        -fMat[1] * invDet,                  // new b = -b/det
        fMat[0] * invDet,                   // new d = a/det
        (fMat[1] * fMat[4] - fMat[0] * fMat[5]) * invDet   // new f = (b*e - a*f)/det
    );
}


// transform points using matrix --
void GMatrix::mapPoints(GPoint dst[], const GPoint src[], int count) const {
    for (int i = 0; i < count; i++) {
        float x = src[i].x;
        float y = src[i].y;
        dst[i].x = fMat[0] * x + fMat[2] * y + fMat[4];
        dst[i].y = fMat[1] * x + fMat[3] * y + fMat[5];
    }
}


