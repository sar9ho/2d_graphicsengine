/*
 *  Copyright 2025 Sarah Ouda
 */

#include "MyEdger.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <optional>

// HELPER: Evaluate a quadratic bezier at param t
static inline GPoint evalQuadPt(const GPoint& a, const GPoint& b, const GPoint& c, float t) {
    float inv = 1 - t;
    return { inv * inv * a.x + 2 * inv * t * b.x + t * t * c.x,
             inv * inv * a.y + 2 * inv * t * b.y + t * t * c.y };
}

// HELPER: same but for cubic bezier
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

// compute # of segments needed to flatten quadratic curve -- uses perpendicular distance from the control point to the chord
static int computeQuadSegments(const GPoint& p0, const GPoint& p1, const GPoint& p2, float tol) {
    float dx = p2.x - p0.x, dy = p2.y - p0.y;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < tol) return 1;
    float cross = std::abs((p1.x - p0.x) * dy - (p1.y - p0.y) * dx);
    float d = cross / len;
    int segs = static_cast<int>(std::ceil(d / tol));
    return std::max(segs, 1);
}

// computes # of segments needed to flatten a cubic curve-- uses the max distance of the middle control points from the chord
static int computeCubicSegments(const GPoint& p0, const GPoint& p1, const GPoint& p2, const GPoint& p3, float tol) {
    float dx = p3.x - p0.x, dy = p3.y - p0.y;
    float len = std::sqrt(dx * dx + dy * dy);
    if (len < tol) return 1;
    float d1 = std::abs((p1.x - p0.x) * dy - (p1.y - p0.y) * dx) / len;
    float d2 = std::abs((p2.x - p0.x) * dy - (p2.y - p0.y) * dx) / len;
    float dmax = std::max(d1, d2);
    int segs = static_cast<int>(std::ceil(dmax / tol)); 
    return std::max(segs, 1);
}

MyEdger::MyEdger(const GPath& path)
    : fEdger(path), fCurrentFlatIndex(0)
{
    // constructed using prof's edger constructor
    // flattened buffer initially empty...
}

std::optional<GPathVerb> MyEdger::next(GPoint pts[]) {
    const float kTolerance = 0.25f;

    // finish emitting flattened segment
    if (!fFlattenedPts.empty() && fCurrentFlatIndex < static_cast<int>(fFlattenedPts.size()) - 1) {
        pts[0] = fFlattenedPts[fCurrentFlatIndex];
        pts[1] = fFlattenedPts[fCurrentFlatIndex + 1];
        fLastPt = pts[1];
        fCurrentFlatIndex++;
        return kLine;
    }

    // if has pending close, emit now
    if (fHasPendingClose) {
        fHasPendingClose = false;
        if (fLastPt != fContourStart) {
            pts[0] = fLastPt;
            pts[1] = fContourStart;
            fLastPt = fContourStart;
            return kLine;
        }
    }

    fFlattenedPts.clear();
    fCurrentFlatIndex = 0;

    std::optional<GPathVerb> optVerb = fEdger.next(pts);
    if (!optVerb.has_value()) {
        return {};
    }

    GPathVerb verb = optVerb.value();

    if (verb == kLine) {
        fLastPt = pts[1];
        return kLine;
    }

    if (verb == kQuad) {
        GPoint p0 = pts[0], p1 = pts[1], p2 = pts[2];
        int segs = computeQuadSegments(p0, p1, p2, kTolerance);
        fFlattenedPts.resize(segs + 1);
        float dt = 1.0f / segs;
        float t = 0.0f;
        for (int i = 0; i <= segs; ++i) {
            fFlattenedPts[i] = evalQuadPt(p0, p1, p2, t);
            t += dt;
        }
        fCurrentFlatIndex = 0;
        pts[0] = fFlattenedPts[0];
        pts[1] = fFlattenedPts[1];
        fLastPt = pts[1];
        fCurrentFlatIndex++;
        return kLine;
    }

    if (verb == kCubic) {
        GPoint p0 = pts[0], p1 = pts[1], p2 = pts[2], p3 = pts[3];
        int segs = computeCubicSegments(p0, p1, p2, p3, kTolerance);
        fFlattenedPts.resize(segs + 1);
        float dt = 1.0f / segs;
        float t = 0.0f;
        for (int i = 0; i <= segs; ++i) {
            fFlattenedPts[i] = evalCubicPt(p0, p1, p2, p3, t);
            t += dt;
        }
        fCurrentFlatIndex = 0;
        pts[0] = fFlattenedPts[0];
        pts[1] = fFlattenedPts[1];
        fLastPt = pts[1];
        fCurrentFlatIndex++;
        return kLine;
    }

    if (verb == kMove) {
        fContourStart = pts[0];
        fLastPt = pts[0];
        fHasPendingClose = true;  // start ofnew contour
        return next(pts);  // skip move, go straight to next edge
    }

    return {};
}
