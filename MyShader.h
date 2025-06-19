/*
 *  Copyright 2025 Sarah Ouda
 */


#include "include/GShader.h"
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPoint.h"
#include "include/GPixel.h"
#include <memory>
#include <vector>

class BitmapShader : public GShader {
public:
    BitmapShader(const GBitmap& bm, const GMatrix& lm);
    BitmapShader(const GBitmap& bm, const GMatrix& lm, GTileMode modeX, GTileMode modeY);

    bool        isOpaque() override;
    bool        setContext(const GMatrix& ctm) override;
    void        shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    GBitmap     bitmap;
    GMatrix     localMatrix;
    GMatrix     inverseMatrix;
    GPoint      deltaX, deltaY;
    bool        opaque;
    GTileMode   fTileModeX, fTileModeY;

    float       applyTileMode(float coord, int max, GTileMode mode);
};

class TriColorShader : public GShader {
public:
    TriColorShader(const GPoint& p0, const GPoint& p1, const GPoint& p2,
                   const GColor& c0, const GColor& c1, const GColor& c2);

    bool        isOpaque() override;
    bool        setContext(const GMatrix& ctm) override;
    void        shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    GPoint      P0, P1, P2;
    GColor      C0, C1, C2;
    float       inv00, inv01, inv10, inv11;
    GPoint      D0;
};

class ProxyShader : public GShader {
public:
    ProxyShader(const GPoint& p0, const GPoint& p1, const GPoint& p2,
                const GPoint& t0, const GPoint& t1, const GPoint& t2,
                GShader* realShader);

    bool        isOpaque() override;
    bool        setContext(const GMatrix& ctm) override;
    void        shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    GShader*    fReal;
    GMatrix     fMap;       
    GPoint      P0, P1, P2;
    GPoint      T0, T1, T2;
};

class CombinedShader : public GShader {
public:
    CombinedShader(std::unique_ptr<GShader> colorShader,
                   std::unique_ptr<GShader> texShader);

    bool        isOpaque() override;
    bool        setContext(const GMatrix& ctm) override;
    void        shadeRow(int x, int y, int count, GPixel row[]) override;

private:
    std::unique_ptr<GShader> fColor, fTex;
    std::vector<GPixel>      fColorBuf, fTexBuf;
};

std::unique_ptr<GShader> makeTriangleShader(
    const GPoint& P0, const GPoint& P1, const GPoint& P2,
    const GColor& C0, const GColor& C1, const GColor& C2,
    const GPoint& T0, const GPoint& T1, const GPoint& T2,
    GShader* baseShader);