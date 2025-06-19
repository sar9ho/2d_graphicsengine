/*
 *  Copyright 2025 Sarah Ouda
 */

#include "starter_canvas.h"
#include "include/GShader.h"
#include "include/GMatrix.h"
#include "include/GPath.h"
#include "MyShader.h" 
#include <algorithm>
#include <vector>
#include <optional>
#include "GPixelUtils.h"
#include "include/GPathBuilder.h"
#include "MyEdger.h"



void MyCanvas::clear(const GColor& color) {

    GPixel fincolor = MultandScale(color);

    for (int y = 0; y < fDevice.height(); ++y) {
        for (int x = 0; x < fDevice.width(); ++x) {
            *fDevice.getAddr(x, y) = fincolor;
        }
    }
}

void MyCanvas::save() {
    fStack.push_back(fCTM); // Save current CTM
}

void MyCanvas::restore() {
    if (!fStack.empty()) {
        fCTM = fStack.back(); // restore last saved CTM
        fStack.pop_back();
    }
}

void MyCanvas::concat(const GMatrix& m) {
    fCTM = GMatrix::Concat(fCTM, m); // multiply current matrix by the new one
}


// help 
static GPixel blendPixel(GPixel src, GPixel dst, GBlendMode blendMode) {
    int srcA = GPixel_GetA(src);
    int srcR = GPixel_GetR(src);
    int srcG = GPixel_GetG(src);
    int srcB = GPixel_GetB(src);

    int dstA = GPixel_GetA(dst);
    int dstR = GPixel_GetR(dst);
    int dstG = GPixel_GetG(dst);
    int dstB = GPixel_GetB(dst);

    int outA, outR, outG, outB;

    switch (blendMode) {
        case GBlendMode::kDst:
            outA = dstA;
            outR = dstR;
            outG = dstG;
            outB = dstB;
            break;

        case GBlendMode::kSrcOver: {
            int invSrcA = 255 - srcA;
            outA = srcA + ((dstA * invSrcA + 127) / 255);
            outR = srcR + ((dstR * invSrcA + 127) / 255);
            outG = srcG + ((dstG * invSrcA + 127) / 255);
            outB = srcB + ((dstB * invSrcA + 127) / 255);
            break;
        }

        case GBlendMode::kSrc:
            outA = srcA;
            outR = srcR;
            outG = srcG;
            outB = srcB;
            break;

        case GBlendMode::kClear:
            // out = 0
            outA = outR = outG = outB = 0;
            break;

        case GBlendMode::kDstOver:
            outA = dstA + ((srcA * (255 - dstA) + 127) / 255);
            outR = dstR + ((srcR * (255 - dstA) + 127) / 255);
            outG = dstG + ((srcG * (255 - dstA) + 127) / 255);
            outB = dstB + ((srcB * (255 - dstA) + 127) / 255);
            break;

        case GBlendMode::kSrcIn:
            outA = (srcA * dstA + 127) / 255;
            outR = (srcR * dstA + 127) / 255;
            outG = (srcG * dstA + 127) / 255;
            outB = (srcB * dstA + 127) / 255;
            break;

        case GBlendMode::kDstIn:
            outA = (dstA * srcA + 127) / 255;
            outR = (dstR * srcA + 127) / 255;
            outG = (dstG * srcA + 127) / 255;
            outB = (dstB * srcA + 127) / 255;
            break;

        case GBlendMode::kSrcOut:
            outA = (srcA * (255 - dstA) + 127) / 255;
            outR = (srcR * (255 - dstA) + 127) / 255;
            outG = (srcG * (255 - dstA) + 127) / 255;
            outB = (srcB * (255 - dstA) + 127) / 255;
            break;

        case GBlendMode::kDstOut:
            outA = (dstA * (255 - srcA) + 127) / 255;
            outR = (dstR * (255 - srcA) + 127) / 255;
            outG = (dstG * (255 - srcA) + 127) / 255;
            outB = (dstB * (255 - srcA) + 127) / 255;
            break;

        case GBlendMode::kSrcATop:
            outA = ((dstA * srcA) + 127) / 255 + ((dstA * (255 - srcA) + 127) / 255);
            outR = ((dstA * srcR) + 127) / 255 + ((dstR * (255 - srcA) + 127) / 255);
            outG = ((dstA * srcG) + 127) / 255 + ((dstG * (255 - srcA) + 127) / 255);
            outB = ((dstA * srcB) + 127) / 255 + ((dstB * (255 - srcA) + 127) / 255);
            break;

        case GBlendMode::kDstATop:
            outA = ((srcA * dstA) + 127) / 255 + ((srcA * (255 - dstA) + 127) / 255);
            outR = ((srcA * dstR) + 127) / 255 + ((srcR * (255 - dstA) + 127) / 255);
            outG = ((srcA * dstG) + 127) / 255 + ((srcG * (255 - dstA) + 127) / 255);
            outB = ((srcA * dstB) + 127) / 255 + ((srcB * (255 - dstA) + 127) / 255);
            break;

        case GBlendMode::kXor:
            outA = (((255 - srcA) * dstA + 127) + ((255 - dstA) * srcA + 127)) / 255;
            outR = (((255 - srcA) * dstR + 127) + ((255 - dstA) * srcR + 127)) / 255;
            outG = (((255 - srcA) * dstG + 127) + ((255 - dstA) * srcG + 127)) / 255;
            outB = (((255 - srcA) * dstB + 127) + ((255 - dstA) * srcB + 127)) / 255;
            break;

        default: {
            int invSrcA = 255 - srcA;
            outA = srcA + ((dstA * invSrcA + 127) / 255);
            outR = srcR + ((dstR * invSrcA + 127) / 255);
            outG = srcG + ((dstG * invSrcA + 127) / 255);
            outB = srcB + ((dstB * invSrcA + 127) / 255);
            break;
        }
    }
    return GPixel_PackARGB(outA, outR, outG, outB);
}



//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------


// Helper to store intersection x and its winding contribution
struct EdgeIntersection {
    float x;
    int winding;
};

void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {
    std::shared_ptr<GPath> tPath = path.transform(fCTM);
    GRect bounds = tPath->bounds();

    // Clamp
    int left   = std::max(0,                GRoundToInt(std::floor(bounds.left)));
    int right  = std::min(fDevice.width(),  GRoundToInt(std::ceil(bounds.right)));
    int top    = std::max(0,                GRoundToInt(std::floor(bounds.top)));
    int bottom = std::min(fDevice.height(), GRoundToInt(std::ceil(bounds.bottom)));

    if (left >= right || top >= bottom) {
        return;
    }

    // extract shader + blend mode
    GBlendMode blendMode = paint.getBlendMode();
    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr);
    if (useShader && !shader->setContext(fCTM)) {
        return;
    }

    // prepare a buffer for shader output
    const int rowWidth = right - left;
    std::vector<GPixel> shaderRow(rowWidth);

    // reserve a vector for intersections to avoid repeated allocations
    std::vector<EdgeIntersection> intersections;
    intersections.reserve(16); // reserve for a typical number of intersections

    // loop over each scanline in bounding box
    for (int y = top; y < bottom; ++y) {
        float scanY = y + 0.5f;  // center of the pixel row
        intersections.clear();

        // create edger to iterate over the path edges
        MyEdger edger(*tPath);
        GPoint pts[GPath::kMaxNextPoints];
        while (std::optional<GPathVerb> verbOpt = edger.next(pts)) {
            // each edge is now a line segment from pts[0] to pts[1]
            GPoint a = pts[0];
            GPoint b = pts[1];

            // determine vertical range of the segment
            float yMin = std::min(a.y, b.y);
            float yMax = std::max(a.y, b.y);

            // if scanline intersects the vertical span of the edge, compute x-intersection
            if (scanY >= yMin && scanY < yMax && (b.y - a.y) != 0) {
                float t = (scanY - a.y) / (b.y - a.y);
                float xIntersect = a.x + t * (b.x - a.x);
                // winding: +1 if the edge goes upward, -1 if downward
                int winding = (b.y > a.y) ? +1 : -1;
                intersections.push_back({ xIntersect, winding });
            }
        }

        // if no intersections on scanline, move on
        if (intersections.empty()) {
            continue;
        }
        
        // sort intersections by x value
        std::sort(intersections.begin(), intersections.end(),
                  [](const EdgeIntersection& a, const EdgeIntersection& b) {
                      return a.x < b.x;
                  });
        
        // use winding rule to determine spans to fill
        int windingAcc = 0;
        int spanStart = -1;
        for (size_t i = 0; i < intersections.size(); ++i) {
            // compute the rounded x coordinate once per intersection
            int xIntersectRounded = GRoundToInt(intersections[i].x);
            windingAcc += intersections[i].winding;
            // start a span when winding becomes non-zero
            if (windingAcc != 0 && spanStart < 0) {
                spanStart = xIntersectRounded;
            }
            // end a span when winding returns to zero
            if (windingAcc == 0 && spanStart >= 0) {
                int spanEnd = xIntersectRounded;
                spanStart = std::max(spanStart, left);
                spanEnd   = std::min(spanEnd, right);
                if (spanStart < spanEnd) {
                    // fill span
                    GPixel* rowPtr = fDevice.getAddr(0, y);
                    if (useShader) {
                        shader->shadeRow(spanStart, y, spanEnd - spanStart, shaderRow.data());
                        bool opaque = shader->isOpaque();
                        if (blendMode == GBlendMode::kSrc || (blendMode == GBlendMode::kSrcOver && opaque)) {
                            for (int x = spanStart; x < spanEnd; ++x) {
                                rowPtr[x] = shaderRow[x - spanStart];
                            }
                        } else {
                            for (int x = spanStart; x < spanEnd; ++x) {
                                GPixel src = shaderRow[x - spanStart];
                                GPixel dst = rowPtr[x];
                                rowPtr[x] = blendPixel(src, dst, blendMode);
                            }
                        }
                    } else {
                        GPixel srcColor = MultandScale(paint.getColor());
                        int srcA = GPixel_GetA(srcColor);
                        bool opaque = (srcA == 255);
                        if (blendMode == GBlendMode::kSrc || (blendMode == GBlendMode::kSrcOver && opaque)) {
                            for (int x = spanStart; x < spanEnd; ++x) {
                                rowPtr[x] = srcColor;
                            }
                        } else {
                            for (int x = spanStart; x < spanEnd; ++x) {
                                GPixel dst = rowPtr[x];
                                rowPtr[x] = blendPixel(srcColor, dst, blendMode);
                            }
                        }
                    }
                }
                // reset span start for subsequent spans
                spanStart = -1;
            }
        }
    }
}



//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------
//PA4------------------------------------------------------------------------------------------

void MyCanvas::drawRect(const GRect& rect, const GPaint& paint) {
    // exit early if rect empty
    if (rect.isEmpty()) {
        return;
    }

    GPoint corners[4] = {
        { rect.left,  rect.top    },
        { rect.right, rect.top    },
        { rect.right, rect.bottom },
        { rect.left,  rect.bottom },
    };
    fCTM.mapPoints(corners, 4); // transform rect's four corners to device space

    // compute device bounding box
    float leftF   = corners[0].x;
    float rightF  = corners[0].x;
    float topF    = corners[0].y;
    float bottomF = corners[0].y;
    for (int i = 1; i < 4; ++i) {
        leftF   = std::min(leftF,   corners[i].x);
        rightF  = std::max(rightF,  corners[i].x);
        topF    = std::min(topF,    corners[i].y);
        bottomF = std::max(bottomF, corners[i].y);
    }

    // float bounds to integer pixel bounds, clamped to the device
    int left   = std::max(0,                GRoundToInt(std::floor(leftF)));
    int right  = std::min(fDevice.width(),  GRoundToInt(std::ceil(rightF)));
    int top    = std::max(0,                GRoundToInt(std::floor(topF)));
    int bottom = std::min(fDevice.height(), GRoundToInt(std::ceil(bottomF)));

    if (left >= right || top >= bottom) {
        return;
    }

    // extract blend mode and shader
    GBlendMode blendMode = paint.getBlendMode();
    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr);

    // opt: if kClear, just zero out the bounding box
    if (blendMode == GBlendMode::kClear) {
        for (int y = top; y < bottom; ++y) {
            GPixel* rowPtr = fDevice.getAddr(0, y);
            for (int x = left; x < right; ++x) {
                rowPtr[x] = 0;  // ARGB = 0
            }
        }
        return;
    }

    // invert the CTM to map device pixels -> local coordinates
    std::optional<GMatrix> maybeInv = fCTM.invert();
    if (!maybeInv.has_value()) {
        return;
    }
    GMatrix inv = maybeInv.value();

    // if we have a shader, set its context
    if (useShader) {
        if (!shader->setContext(fCTM)) {
            return;
        }
    }

    // local-space bounds for inside checks:
    const float l = rect.left;
    const float r = rect.right;
    const float t = rect.top;
    const float b = rect.bottom;
    
    const int rowWidth = right - left;
    std::vector<GPixel> row(rowWidth);    // buffer for shader output if used

    // loop over each row in the bounding box
    for (int y = top; y < bottom; ++y) {
        if (useShader) {
            shader->shadeRow(left, y, rowWidth, row.data());
        }

        GPixel* deviceRow = fDevice.getAddr(0, y);

        // --- Forward Differencing for local coordinate computation ---------------
        // compute local coordinate for leftmost pixel center
        GPoint localStart;
        {
            GPoint devicePoint = { left + 0.5f, y + 0.5f };
            inv.mapPoints(&localStart, &devicePoint, 1);
        }
        // compute local coordinate for the next pixel center
        GPoint next;
        {
            GPoint devicePoint = { left + 1 + 0.5f, y + 0.5f };
            inv.mapPoints(&next, &devicePoint, 1);
        }
        // per-pixel delta in local space
        GPoint localDelta = { next.x - localStart.x, next.y - localStart.y };

        // use cumulative addition to compute each pixel's local coordinate
        GPoint local = localStart;
        for (int x = left; x < right; ++x) {
            // check if the computed local coordinate is inside the original rect
            if (local.x >= l && local.x < r && local.y >= t && local.y < b) {
                GPixel srcPixel;
                if (!useShader) {
                    // colid color case
                    srcPixel = MultandScale(paint.getColor());
                } else {
                    // use shader's precomputed row
                    srcPixel = row[x - left];
                }
                // blend into destination
                GPixel dstPixel = deviceRow[x];
                deviceRow[x] = blendPixel(srcPixel, dstPixel, blendMode);
            }
            // increment local coordinate for the next pixel
            local.x += localDelta.x;
            local.y += localDelta.y;
        }
    }
}


void MyCanvas::drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) {
    // polyg must have at least 3 points
    if (count < 3) return;

    // transform polygon points using the current CTM
    std::vector<GPoint> transformedPoints(count);
    fCTM.mapPoints(transformedPoints.data(), points, count);

    // compute bounding box of the transformed polygon
    GRect bounds = { transformedPoints[0].x, transformedPoints[0].y,
                     transformedPoints[0].x, transformedPoints[0].y };
    for (int i = 1; i < count; i++) {
        bounds.left   = std::min(bounds.left, transformedPoints[i].x);
        bounds.right  = std::max(bounds.right, transformedPoints[i].x);
        bounds.top    = std::min(bounds.top, transformedPoints[i].y);
        bounds.bottom = std::max(bounds.bottom, transformedPoints[i].y);
    }

    // clamp bounding box to canvas
    int top    = std::max(GRoundToInt(bounds.top), 0);
    int bottom = std::min(GRoundToInt(bounds.bottom), fDevice.height());

    // reserve space for intersection x-values
    std::vector<int> xIntersections;
    xIntersections.reserve(count);

    // retrieve the shader from the paint
    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr);
    //  get the blend mode
    GBlendMode blendMode = paint.getBlendMode();

    // If shader attached, prepare its context
    if (useShader) {
        if (!shader->setContext(fCTM)) {
            return;
        }
    }

    // iterate over each scanline within bounding box
    for (int y = top; y < bottom; y++) {
        xIntersections.clear();
        float scanlineY = y + 0.5f;

        // compute the x-intersections with polygon edges
        for (int i = 0; i < count; i++) {
            int j = (i + 1) % count;
            float y1 = transformedPoints[i].y, y2 = transformedPoints[j].y;
            float x1 = transformedPoints[i].x, x2 = transformedPoints[j].x;

            if ((y1 <= scanlineY && y2 > scanlineY) ||
                (y2 <= scanlineY && y1 > scanlineY)) {
                float t = (scanlineY - y1) / (y2 - y1);
                float xIntersect = x1 + t * (x2 - x1);
                xIntersections.push_back(GRoundToInt(xIntersect));
            }
        }
        std::sort(xIntersections.begin(), xIntersections.end());

        // get pointer to current scanline
        GPixel* rowPtr = fDevice.getAddr(0, y);

        // process spans btwn intersection pairs
        for (size_t i = 0; i + 1 < xIntersections.size(); i += 2) {
            int xStart = std::max(xIntersections[i], 0);
            int xEnd   = std::min(xIntersections[i + 1], fDevice.width());

            if (useShader) {
                // use shader to generate pixel colors for the span
                std::vector<GPixel> shaderPixels(xEnd - xStart);
                shader->shadeRow(xStart, y, xEnd - xStart, shaderPixels.data());
                bool shaderOpaque = shader->isOpaque();
                // if blend mode is kSrc or kSrcOver + shader output opaque, write shader output directly
                if (blendMode == GBlendMode::kSrc || (blendMode == GBlendMode::kSrcOver && shaderOpaque)) {
                    for (int x = xStart; x < xEnd; x++) {
                        rowPtr[x] = shaderPixels[x - xStart];
                    }
                } else {
                    // else, blend each shader pixel w destination
                    for (int x = xStart; x < xEnd; x++) {
                        GPixel srcPixel = shaderPixels[x - xStart];
                        GPixel dst = rowPtr[x];
                        rowPtr[x] = blendPixel(srcPixel, dst, blendMode);
                    }
                }
            } else {
                // no shader: use solid color from paint
                GPixel srcColor = MultandScale(paint.getColor());
                int srcA = GPixel_GetA(srcColor);
                bool opaque = (srcA == 255);
                // if blend mode is kSrc or kSrcOver with opaque source, write directly
                if (blendMode == GBlendMode::kSrc || (blendMode == GBlendMode::kSrcOver && opaque)) {
                    for (int x = xStart; x < xEnd; x++) {
                        rowPtr[x] = srcColor;
                    }
                } else {
                    // otherwise, blend per pixel using blendPixel
                    for (int x = xStart; x < xEnd; x++) {
                        GPixel dst = rowPtr[x];
                        rowPtr[x] = blendPixel(srcColor, dst, blendMode);
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------------

inline GPoint bilerpPoint(const GPoint verts[4], float u, float v) {
    return (1 - u) * (1 - v) * verts[0] + u * (1 - v) * verts[1] + (1 - u) * v * verts[3] +
           u * v * verts[2];
}

inline GColor bilerpColor(const GColor colors[4], float u, float v) {
    return (1 - u) * (1 - v) * colors[0] + u * (1 - v) * colors[1] + (1 - u) * v * colors[3] +
           u * v * colors[2];
}


void MyCanvas::drawMesh(
    const GPoint verts[],
    const GColor colors[],
    const GPoint texs[],
    int count,
    const int indices[],
    const GPaint& paint) {

    bool haveColors = colors != nullptr;
    bool haveTexs   = texs   != nullptr;

    for (int t = 0; t < count; ++t) {
        int n = 3 * t;
        int i0 = indices[n + 0], i1 = indices[n + 1], i2 = indices[n + 2];

        GPoint P0 = verts[i0], P1 = verts[i1], P2 = verts[i2];
        GColor C0 = haveColors ? colors[i0] : GColor{0,0,0,0},
               C1 = haveColors ? colors[i1] : GColor{0,0,0,0},
               C2 = haveColors ? colors[i2] : GColor{0,0,0,0};
        GPoint T0 = haveTexs ? texs[i0] : GPoint{0,0},
               T1 = haveTexs ? texs[i1] : GPoint{0,0},
               T2 = haveTexs ? texs[i2] : GPoint{0,0};

        // Winding check
        float dx1 = P1.x - P0.x, dy1 = P1.y - P0.y;
        float dx2 = P2.x - P0.x, dy2 = P2.y - P0.y;
        if (dx1*dy2 - dy1*dx2 < 0) {
            std::swap(P1, P2);
            std::swap(C1, C2);
            std::swap(T1, T2);
        }

        GPaint triPaint = paint;

        if (haveColors && haveTexs) {
            auto colorShader = std::make_unique<TriColorShader>(P0, P1, P2, C0, C1, C2);
            // GPoint triPts[3] = { P0, P1, P2 };
            // GPoint triTexs[3] = { T0, T1, T2 };
            auto proxyShader = std::make_unique<ProxyShader>(P0, P1, P2, T0, T1, T2, paint.peekShader());
            auto composedShader = std::make_unique<CombinedShader>(std::move(colorShader), std::move(proxyShader));
            triPaint.setShader(std::move(composedShader));

        } else if (haveColors) {
            auto colorShader = std::make_unique<TriColorShader>(P0, P1, P2, C0, C1, C2);
            triPaint.setShader(std::move(colorShader));
        } else if (haveTexs) {
            auto proxyShader = std::make_unique<ProxyShader>(P0, P1, P2, T0, T1, T2, paint.peekShader());
            triPaint.setShader(std::move(proxyShader));
        } else {
            // no shader --> just plain paint color
        }

        GPoint triPts[3] = { P0, P1, P2 };
        this->drawConvexPolygon(triPts, 3, triPaint);
    }
}


void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                        int level, const GPaint& paint) {
    GPoint subVerts[4];
    GColor subColors[4];
    GPoint subTexs[4];


    int indices[6] = {0, 1, 3, 1, 2, 3};
    for (int u = 0; u < level + 1; ++u) {
        auto u0 = static_cast<float>(u) / (level + 1);
        auto u1 = static_cast<float>(u + 1) / (level + 1);
        for (int v = 0; v < level + 1; ++v) {
            auto v0 = static_cast<float>(v) / (level + 1);
            auto v1 = static_cast<float>(v + 1) / (level + 1);

            float uv[4][2] = {{u0, v0}, {u1, v0}, {u1, v1}, {u0, v1}};

            for (int i = 0; i < 4; ++i) {
                subVerts[i] = bilerpPoint(verts, uv[i][0], uv[i][1]);
                if (colors != nullptr) {
                    subColors[i] = bilerpColor(colors, uv[i][0], uv[i][1]);
                }
                if (texs != nullptr) {
                    subTexs[i] = bilerpPoint(texs, uv[i][0], uv[i][1]);
                }
            }

            drawMesh(subVerts, colors != nullptr ? subColors : nullptr,
                     texs != nullptr ? subTexs : nullptr, 2, indices, paint);
        }
    }
}



std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& device) {
    return std::unique_ptr<GCanvas>(new MyCanvas(device));
}


std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    canvas->clear({1, 1, 1, 1});

    GPathBuilder pb;
    GPoint center = { dim.width / 2.0f, dim.height / 2.0f };
    float radius = std::min(dim.width, dim.height) * 0.3f;

    pb.addCircle(center, radius, GPathDirection::kCW);

    GPaint paint;
    paint.setColor({1, 0, 0, 1});  // red

    canvas->drawPath(*pb.detach(), paint);  

    return "Draw red circle using addCircle";
}
