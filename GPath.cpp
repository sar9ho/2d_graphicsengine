/*
 *  Copyright 2025 Sarah Ouda
 */

 #include "include/GPath.h"
 #include "include/GMatrix.h"
 #include <algorithm>
 #include <optional>
 #include <cfloat>
 #include <cmath>
 
 // Inline helper: Evaluate a quadratic Bezier at parameter t
 inline float evalQuad(float p0, float p1, float p2, float t) {
     float inv = 1 - t;
     return inv * inv * p0 + 2 * inv * t * p1 + t * t * p2;
 }
 
 // Inline helper: Evaluate a cubic Bezier at parameter t
 inline float evalCubic(float p0, float p1, float p2, float p3, float t) {
     float inv = 1 - t;
     return inv * inv * inv * p0 +
            3 * inv * inv * t * p1 +
            3 * inv * t * t * p2 +
            t * t * t * p3;
 }
 
 GRect GPath::bounds() const {
     if (fVbs.empty()) {
         return GRect{0, 0, 0, 0};
     }
 
     float minX = FLT_MAX, minY = FLT_MAX;
     float maxX = -FLT_MAX, maxY = -FLT_MAX;
     size_t ptIndex = 0;
     GPoint current; // current starting point for the segment
 
     // Iterate through each segment in the path using fVbs
     for (auto verb : fVbs) {
         switch (verb) {
             case kMove: {
                 current = fPts[ptIndex++];
                 minX = std::min(minX, current.x);
                 maxX = std::max(maxX, current.x);
                 minY = std::min(minY, current.y);
                 maxY = std::max(maxY, current.y);
                 break;
             }
             case kLine: {
                 GPoint p = fPts[ptIndex++];
                 minX = std::min({minX, current.x, p.x});
                 maxX = std::max({maxX, current.x, p.x});
                 minY = std::min({minY, current.y, p.y});
                 maxY = std::max({maxY, current.y, p.y});
                 current = p;
                 break;
             }
             case kQuad: {
                 // Quadratic: three points (p0 = current, p1, p2)
                 GPoint p1 = fPts[ptIndex++];
                 GPoint p2 = fPts[ptIndex++];
                 minX = std::min({minX, current.x, p2.x});
                 maxX = std::max({maxX, current.x, p2.x});
                 minY = std::min({minY, current.y, p2.y});
                 maxY = std::max({maxY, current.y, p2.y});
 
                 // For each coordinate, solve for t where the derivative is zero
                 float denomX = current.x - 2 * p1.x + p2.x;
                 if (std::fabs(denomX) > 1e-6f) {
                     float t = (current.x - p1.x) / denomX;
                     if (t > 0 && t < 1) {
                         float x = evalQuad(current.x, p1.x, p2.x, t);
                         minX = std::min(minX, x);
                         maxX = std::max(maxX, x);
                     }
                 }
                 float denomY = current.y - 2 * p1.y + p2.y;
                 if (std::fabs(denomY) > 1e-6f) {
                     float t = (current.y - p1.y) / denomY;
                     if (t > 0 && t < 1) {
                         float y = evalQuad(current.y, p1.y, p2.y, t);
                         minY = std::min(minY, y);
                         maxY = std::max(maxY, y);
                     }
                 }
                 current = p2;
                 break;
             }
             case kCubic: {
                 // Cubic: four points (p0 = current, p1, p2, p3)
                 GPoint p1 = fPts[ptIndex++];
                 GPoint p2 = fPts[ptIndex++];
                 GPoint p3 = fPts[ptIndex++];
                 minX = std::min({minX, current.x, p3.x});
                 maxX = std::max({maxX, current.x, p3.x});
                 minY = std::min({minY, current.y, p3.y});
                 maxY = std::max({maxY, current.y, p3.y});
 
                 // For cubic, the derivative is quadratic.
                 // X derivative: solve for t in: ax * t^2 + bx * t + cx = 0
                 float ax = -current.x + 3 * p1.x - 3 * p2.x + p3.x;
                 float bx = 2 * (current.x - 2 * p1.x + p2.x);
                 float cx = -current.x + p1.x;
                 float rootsX[2];
                 int numRootsX = 0;
                 if (std::fabs(ax) < 1e-6f) {
                     if (std::fabs(bx) > 1e-6f) {
                         rootsX[numRootsX++] = -cx / bx;
                     }
                 } else {
                     float disc = bx * bx - 4 * ax * cx;
                     if (disc >= 0) {
                         float sqrtDisc = std::sqrt(disc);
                         rootsX[numRootsX++] = (-bx + sqrtDisc) / (2 * ax);
                         rootsX[numRootsX++] = (-bx - sqrtDisc) / (2 * ax);
                     }
                 }
                 for (int i = 0; i < numRootsX; ++i) {
                     float t = rootsX[i];
                     if (t > 0 && t < 1) {
                         float x = evalCubic(current.x, p1.x, p2.x, p3.x, t);
                         minX = std::min(minX, x);
                         maxX = std::max(maxX, x);
                     }
                 }
 
                 // Repeat for y:
                 float ay = -current.y + 3 * p1.y - 3 * p2.y + p3.y;
                 float by = 2 * (current.y - 2 * p1.y + p2.y);
                 float cy = -current.y + p1.y;
                 float rootsY[2];
                 int numRootsY = 0;
                 if (std::fabs(ay) < 1e-6f) {
                     if (std::fabs(by) > 1e-6f) {
                         rootsY[numRootsY++] = -cy / by;
                     }
                 } else {
                     float disc = by * by - 4 * ay * cy;
                     if (disc >= 0) {
                         float sqrtDisc = std::sqrt(disc);
                         rootsY[numRootsY++] = (-by + sqrtDisc) / (2 * ay);
                         rootsY[numRootsY++] = (-by - sqrtDisc) / (2 * ay);
                     }
                 }
                 for (int i = 0; i < numRootsY; ++i) {
                     float t = rootsY[i];
                     if (t > 0 && t < 1) {
                         float y = evalCubic(current.y, p1.y, p2.y, p3.y, t);
                         minY = std::min(minY, y);
                         maxY = std::max(maxY, y);
                     }
                 }
                 current = p3;
                 break;
             }
         } // end switch
     } // end for
 
     return GRect{minX, minY, maxX, maxY};
 }
  
 void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
     // Compute the intermediate points.
     GPoint ab = {
         (1 - t) * src[0].x + t * src[1].x,
         (1 - t) * src[0].y + t * src[1].y
     };
     GPoint bc = {
         (1 - t) * src[1].x + t * src[2].x,
         (1 - t) * src[1].y + t * src[2].y
     };
     GPoint abc = {
         (1 - t) * ab.x + t * bc.x,
         (1 - t) * ab.y + t * bc.y
     };
 
     // Store quadratics
     dst[0] = src[0];
     dst[1] = ab;
     dst[2] = abc;  
     dst[3] = bc;
     dst[4] = src[2];
 }
 
 void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
     // first level of lerp
     GPoint ab = {
         (1 - t) * src[0].x + t * src[1].x,
         (1 - t) * src[0].y + t * src[1].y
     };
     GPoint bc = {
         (1 - t) * src[1].x + t * src[2].x,
         (1 - t) * src[1].y + t * src[2].y
     };
     GPoint cd = {
         (1 - t) * src[2].x + t * src[3].x,
         (1 - t) * src[2].y + t * src[3].y
     };
 
     //second level
     GPoint abc = {
         (1 - t) * ab.x + t * bc.x,
         (1 - t) * ab.y + t * bc.y
     };
     GPoint bcd = {
         (1 - t) * bc.x + t * cd.x,
         (1 - t) * bc.y + t * cd.y
     };
 
     // final lerp
     GPoint abcd = {
         (1 - t) * abc.x + t * bcd.x,
         (1 - t) * abc.y + t * bcd.y
     };
 
     // Store subdivided cubics
     dst[0] = src[0];
     dst[1] = ab;
     dst[2] = abc;
     dst[3] = abcd;  // Shared endpoint
     dst[4] = bcd;
     dst[5] = cd;
     dst[6] = src[3];
 }
 
 // Helper functions to evaluate points along curves
 static inline GPoint evalQuadPt(const GPoint& a, const GPoint& b, const GPoint& c, float t) {
     float inv = 1 - t;
     return { inv * inv * a.x + 2 * inv * t * b.x + t * t * c.x,
              inv * inv * a.y + 2 * inv * t * b.y + t * t * c.y };
 }
 
 static inline GPoint evalCubicPt(const GPoint& a, const GPoint& b, const GPoint& c, const GPoint& d, float t) {
     float inv = 1 - t;
     return { inv * inv * inv * a.x +
              3 * inv * inv * t * b.x +
              3 * inv * t * t * c.x +
              t * t * t * d.x,
              inv * inv * inv * a.y +
              3 * inv * inv * t * b.y +
              3 * inv * t * t * c.y +
              t * t * t * d.y };
 }
 