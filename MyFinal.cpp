/*
 *  Copyright 2025 Sarah Ouda
 */

#include "include/GFinal.h"
#include "include/GPoint.h"
#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include "starter_canvas.h"
#include "GPixelUtils.h"
#include <memory>


class MyFinal : public GFinal {
public:

    //HELPERS FOR DRAW QUADRATIC COONS==========================================================
    
    inline GPoint evalQuad(const GPoint& p0, const GPoint& p1, const GPoint& p2, float t) {
        float inv = 1 - t;
        return {
            inv * inv * p0.x + 2 * inv * t * p1.x + t * t * p2.x,
            inv * inv * p0.y + 2 * inv * t * p1.y + t * t * p2.y
        };
    }
    
    inline GPoint lerp(const GPoint& p0, const GPoint& p1, float t) {
        return { (1 - t) * p0.x + t * p1.x,
                 (1 - t) * p0.y + t * p1.y };
    }
    
    inline GPoint bilerpPoint(const GPoint verts[4], float u, float v) {
        return (1 - u) * (1 - v) * verts[0]
             +      u  * (1 - v) * verts[1]
             + (1 - u) *      v  * verts[3]
             +      u  *      v  * verts[2];
    }
    // ---------------------------------------------------------------------------------------

    // Voronoi Shader ========================================================
    // Voronoi Shader ========================================================
    // Voronoi Shader ========================================================

    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[],
                                                 const GColor colors[],
                                                 int count) override {
            class VoronoiShader: public GShader {
            private:
                std::vector<GPoint> myPoints;
                std::vector<GColor> myColors;
                int myCount;
                GMatrix invTransform;
                GMatrix M;
            public:
                VoronoiShader(const GPoint points[], const GColor colors[], int count, GMatrix m){
                    myCount = count;
                    M = m;
                    for(int i = 0; i < count; i++) {
                        myPoints.push_back(points[i]);
                        myColors.push_back(colors[i]);
                    }
                };
                bool isOpaque() override {
                    return false;
                };
                
                bool setContext(const GMatrix& ctm) override {
                    auto inv = (ctm * M).invert();
                    if (inv.has_value()) {
                        invTransform = inv.value();
                        return true;
                    }
                    return false;
                }
                
                void shadeRow(int x, int y, int count, GPixel row[]) override {
                    GPoint p_t = invTransform * GPoint{x + 0.5f, y + 0.5f};
                    float px = p_t.x;
                    float dpx = invTransform.e0().x;


                    for (int i = 0; i < count; i++, px += dpx) {
                        float minDist = std::numeric_limits<float>::max();
                        int minIndex = 0;
                        for (int j = 0; j < myCount; j++) {
                            float dist = (myPoints[j] - GPoint{px, p_t.y}).length();
                            if (dist < minDist) {
                                minDist = dist;
                                minIndex = j;
                            }
                        }
                        row[i] = MultandScale(myColors[minIndex]);
                    }
                };
            };
            GMatrix m = {1, 0, 0, 0, 1, 0 };  
            return std::make_shared<VoronoiShader>(points, colors, count, m);
        }
    // ----------------------------------------------------------------------------------------

    // Stroke Polygon ========================================================================
    // Stroke Polygon ========================================================================
    // Stroke Polygon ========================================================================

    std::shared_ptr<GPath> strokePolygon(const GPoint pts[], int count,
                                     float width, bool isClosed) override {
        if (count < 2 || width <= 0) {
            return nullptr;
        }

        float radius = width * 0.5f;
        auto path = std::make_shared<GPathBuilder>();

        auto unitPerp = [](GPoint p0, GPoint p1) -> GPoint {
            GPoint dir = { p1.x - p0.x, p1.y - p0.y };
            float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
            if (len == 0) return {0,0};
            return { -(dir.y / len), (dir.x / len) };  // 90 degree ccw
        };

        int N = count;
        if (!isClosed) N -= 1;

        for (int i = 0; i < N; ++i) {
            int j = (i + 1) % count;

            GPoint p0 = pts[i];
            GPoint p1 = pts[j];

            GPoint n = unitPerp(p0, p1);
            GPoint n0 = { n.x * radius, n.y * radius };
            GPoint n1 = { -n0.x, -n0.y };

            GPoint v0 = { p0.x + n0.x, p0.y + n0.y };
            GPoint v1 = { p1.x + n0.x, p1.y + n0.y };
            GPoint v2 = { p1.x + n1.x, p1.y + n1.y };
            GPoint v3 = { p0.x + n1.x, p0.y + n1.y };

            path->moveTo(v0);
            path->lineTo(v1);
            path->lineTo(v2);
            path->lineTo(v3);
            // let drawPath handle
        }

        // add round caps for open strokes
        if (!isClosed) {
            GPoint start = pts[0];
            GPoint end = pts[count - 1];
            path->addCircle(start, radius, GPathDirection::kCW);
            path->addCircle(end, radius, GPathDirection::kCW);
        }

        // add round joins at all points
        for (int i = 0; i < count; ++i) {
            path->addCircle(pts[i], radius, GPathDirection::kCW);
        }

        return path->detach();
    }


//------------------------------------------------------------------------------------------------------

// drawQuadraticCoons=================================================================================
// drawQuadraticCoons=================================================================================
// drawQuadraticCoons=================================================================================

    void drawQuadraticCoons(GCanvas* canvas,
        const GPoint pts[8],
        const GPoint tex[4],
        int level,
        const GPaint& paint) {
        int subdiv = level + 1;
        int numVertsPerSide = subdiv + 1;

        std::vector<GPoint> verts;
        std::vector<GPoint> texs;
        std::vector<int> indices;

        // generate grid of points
        for (int j = 0; j < numVertsPerSide; ++j) {
        float v = float(j) / subdiv;
        for (int i = 0; i < numVertsPerSide; ++i) {
        float u = float(i) / subdiv;

        GPoint TB = lerp(evalQuad(pts[0], pts[1], pts[2], u),
            evalQuad(pts[6], pts[5], pts[4], u),
            v);

        GPoint LR = lerp(evalQuad(pts[0], pts[7], pts[6], v),
            evalQuad(pts[2], pts[3], pts[4], v),
            u);

        GPoint corners[4] = { pts[0], pts[2], pts[4], pts[6] };
        GPoint CornerInterp = bilerpPoint(corners, u, v);

        verts.push_back({ TB.x + LR.x - CornerInterp.x, TB.y + LR.y - CornerInterp.y });

        if (tex != nullptr) {
        texs.push_back(bilerpPoint(tex, u, v));
        }
        }
        }

        // create triangles from grid
        for (int j = 0; j < subdiv; ++j) {
        for (int i = 0; i < subdiv; ++i) {
        int idx0 = j * numVertsPerSide + i;
        int idx1 = idx0 + 1;
        int idx2 = idx0 + numVertsPerSide;
        int idx3 = idx2 + 1;

        // 2 triangles r (idx0, idx1, idx2) & (idx1, idx3, idx2)
        indices.push_back(idx0);
        indices.push_back(idx1);
        indices.push_back(idx2);

        indices.push_back(idx1);
        indices.push_back(idx3);
        indices.push_back(idx2);
        }
        }

        // call drawMesh to draw everything
        canvas->drawMesh(
            verts.data(),
            nullptr,
            (tex == nullptr) ? nullptr : texs.data(),
            indices.size() / 3,   // triangle count
            indices.data(),       // array of actual indices
            paint
        );
        }

};


std::unique_ptr<GFinal> GCreateFinal() {
    return std::make_unique<MyFinal>();
}
