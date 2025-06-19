/*
 *  Copyright 2025 Sarah Ouda
 */

 #ifndef _g_starter_canvas_h_
 #define _g_starter_canvas_h_
 
 #include "include/GCanvas.h"
 #include "include/GRect.h"
 #include "include/GColor.h"
 #include "include/GBitmap.h"
 
 #include "include/GMatrix.h"
 #include <vector> 

class MyCanvas : public GCanvas {
public:
    MyCanvas(const GBitmap& device) : fDevice(device), fCTM() {}

    void clear(const GColor& color) override;
    void drawRect(const GRect& rect, const GPaint& paint) override;
    void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override;
    virtual void drawPath(const GPath& path, const GPaint& paint) override;
    
    void drawMesh(const GPoint verts[],
        const GColor colors[],
        const GPoint texs[],
        int count,
        const int indices[],
        const GPaint& paint) override;

    void drawQuad(const GPoint verts[4],
            const GColor colors[4],
            const GPoint texs[4],
            int level,
            const GPaint& paint) override;


    void save() override;
    void restore() override;
    void concat(const GMatrix&) override;

private:
    const GBitmap fDevice;
    GMatrix fCTM;                // current transformation matrix
    std::vector<GMatrix> fStack; // stack for saving/restoring matrices

    // PA6 helper -- rasterize one triangle with per‑vertex color & tex
    void drawTriangle(const GPoint  p0, const GPoint  p1, const GPoint  p2,
        const GColor  c0, const GColor  c1, const GColor  c2,
        const GPoint  t0, const GPoint  t1, const GPoint  t2,
        const GPaint& paint);
};

 
 #endif
 