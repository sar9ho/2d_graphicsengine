/*
 *  Copyright 2025 Sarah Ouda
 */


#include "MyShader.h"
#include <algorithm>
#include "include/GPoint.h"
#include "include/GShader.h" 
#include "include/GBitmap.h"
#include "include/GMatrix.h"
#include "include/GPixel.h"
#include <cmath>


BitmapShader::BitmapShader(const GBitmap& bm, const GMatrix& lm)
    : bitmap(bm), localMatrix(lm), opaque(bm.isOpaque()),
      fTileModeX(GTileMode::kClamp), fTileModeY(GTileMode::kClamp) { }

BitmapShader::BitmapShader(const GBitmap& bm, const GMatrix& lm, GTileMode modeX, GTileMode modeY)
    : bitmap(bm), localMatrix(lm), opaque(bm.isOpaque()),
      fTileModeX(modeX), fTileModeY(modeY) {}

bool BitmapShader::isOpaque() {
    return opaque;
}

bool BitmapShader::setContext(const GMatrix& ctm) {
    GMatrix combined = GMatrix::Concat(ctm, localMatrix);
    auto inv = combined.invert();
    if (!inv) return false;
    inverseMatrix = *inv;
    GPoint origin = inverseMatrix * GPoint{0.5f, 0.5f};
    GPoint dx = inverseMatrix * GPoint{1.5f, 0.5f};
    deltaX = { dx.x - origin.x, dx.y - origin.y };
    return true;
}

void BitmapShader::shadeRow(int x, int y, int count, GPixel row[]) {
    int bmpW = bitmap.width();
    int bmpH = bitmap.height();
    GPoint src = inverseMatrix * GPoint{x + 0.5f, y + 0.5f};
    for (int i = 0; i < count; ++i) {
        float u = applyTileMode(src.x, bmpW, fTileModeX);
        float v = applyTileMode(src.y, bmpH, fTileModeY);
        int sx = std::clamp(int(std::floor(u)), 0, bmpW - 1);
        int sy = std::clamp(int(std::floor(v)), 0, bmpH - 1);
        row[i] = *bitmap.getAddr(sx, sy);
        src.x += deltaX.x;
        src.y += deltaX.y;
    }
}

inline float BitmapShader::applyTileMode(float coord, int max, GTileMode mode) {
    switch (mode) {
        case GTileMode::kClamp:
            return std::clamp(coord, 0.0f, float(max) - 1);
        case GTileMode::kRepeat:
            return float(std::fmod(std::fmod(coord, max) + max, max));
        case GTileMode::kMirror: {
            float m = float(std::fmod(coord, max * 2));
            if (m < 0) m += max * 2;
            return m < max ? m : max * 2 - m;
        }
    }
    return coord;
}


std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& matrix, GTileMode modeX, GTileMode modeY) {
    return std::make_shared<BitmapShader>(bitmap, matrix, modeX, modeY);
}
std::shared_ptr<GShader> GCreateBitmapShader(const GBitmap& bitmap, const GMatrix& matrix, GTileMode mode) {
    return std::make_shared<BitmapShader>(bitmap, matrix, mode, mode);
}

// --- TriColorShader Implementation ---
// --- TriColorShader Implementation ---

TriColorShader::TriColorShader(const GPoint& p0, const GPoint& p1, const GPoint& p2,
    const GColor& c0, const GColor& c1, const GColor& c2)
: P0(p0), P1(p1), P2(p2), C0(c0), C1(c1), C2(c2) {
// no longer need to compute inv00…inv11 here,,, --> will be done in setContext()
}

bool TriColorShader::isOpaque() {
    return C0.a == 1 && C1.a == 1 && C2.a == 1;
}

bool TriColorShader::setContext(const GMatrix& ctm) {
    // transform triangle vertices into device space
    GPoint d0 = P0, d1 = P1, d2 = P2;
    ctm.mapPoints(&d0, 1);
    ctm.mapPoints(&d1, 1);
    ctm.mapPoints(&d2, 1);

    // recompute the inverse of the 2×2 matrix [ d1 - d0 | d2 - d0 ]
    float v1x = d1.x - d0.x, v1y = d1.y - d0.y;
    float v2x = d2.x - d0.x, v2y = d2.y - d0.y;
    float det  = v1x * v2y - v2x * v1y;
    if (det == 0) return false;  // if triangle degenerate
    float invDet = 1.0f / det;

    inv00 =  invDet * v2y;
    inv01 = -invDet * v2x;
    inv10 = -invDet * v1y;
    inv11 =  invDet * v1x;

    D0 = d0; 

    return true;
}

void TriColorShader::shadeRow(int x, int y, int count, GPixel row[]) {
    // premul vertex colors
    float r0 = C0.r * C0.a,  g0 = C0.g * C0.a,  b0 = C0.b * C0.a,  a0 = C0.a;
    float r1 = C1.r * C1.a,  g1 = C1.g * C1.a,  b1 = C1.b * C1.a,  a1 = C1.a;
    float r2 = C2.r * C2.a,  g2 = C2.g * C2.a,  b2 = C2.b * C2.a,  a2 = C2.a;

    for (int i = 0; i < count; ++i) {
        float sx = x + i + 0.5f,
        sy = y     + 0.5f;

        // computation of  barycentrics relative to device‐space D0
        float dx    = sx - D0.x,
        dy    = sy - D0.y;
        float beta  = inv00 * dx + inv01 * dy;
        float gamma = inv10 * dx + inv11 * dy;
        float alpha = 1 - beta - gamma;

        // interpolate premul
        float a = alpha*a0 + beta*a1 + gamma*a2;
        float r = alpha*r0 + beta*r1 + gamma*r2;
        float g = alpha*g0 + beta*g1 + gamma*g2;
        float b = alpha*b0 + beta*b1 + gamma*b2;

        // clamp
        a = std::clamp(a, 0.0f, 1.0f);
        r = std::clamp(r, 0.0f, a);
        g = std::clamp(g, 0.0f, a);
        b = std::clamp(b, 0.0f, a);

        // ..to 8‑bit
        int A = GRoundToInt(a * 255),
        R = GRoundToInt(r * 255),
        G = GRoundToInt(g * 255),
        B = GRoundToInt(b * 255);
        R = std::min(R, A);
        G = std::min(G, A);
        B = std::min(B, A);

        row[i] = GPixel_PackARGB(A, R, G, B); //pack
    }
}



// --- ProxyShader Implementation ---

ProxyShader::ProxyShader(const GPoint& p0, const GPoint& p1, const GPoint& p2,
    const GPoint& t0, const GPoint& t1, const GPoint& t2,
    GShader* realShader)
: fReal(realShader)
, P0(p0),  P1(p1),  P2(p2)
, T0(t0),  T1(t1),  T2(t2)
{}

// j forward opacity
bool ProxyShader::isOpaque() {
    return fReal->isOpaque();
}


bool ProxyShader::setContext(const GMatrix& ctm) {
    // bring triangle into d space
    GPoint d0 = P0, d1 = P1, d2 = P2;
    ctm.mapPoints(&d0, 1);
    ctm.mapPoints(&d1, 1);
    ctm.mapPoints(&d2, 1);

    // build barycentric→device matrix (triDev)
    GVector dv0 = { d1.x - d0.x, d1.y - d0.y };
    GVector dv1 = { d2.x - d0.x, d2.y - d0.y };
    GMatrix triDev(dv0, dv1, d0);

    // build texture→barycentric inverse (T_inv) 
    float u1x = T1.x - T0.x,  u1y = T1.y - T0.y,
          u2x = T2.x - T0.x,  u2y = T2.y - T0.y;

    float detT   = u1x*u2y - u2x*u1y;
    if (detT == 0) return false;
    float invDet = 1.0f / detT;

    // 2×2 inverse M = 1/det * [ u2y  -u2x ;  -u1y  u1x ]
    float inv00 =  invDet * u2y,
          inv01 = -invDet * u2x,
          inv10 = -invDet * u1y,
          inv11 =  invDet * u1x;

    // build gmatrix from columns of M
    GVector col0{ inv00, inv10 };   // first column
    GVector col1{ inv01, inv11 };   // second column

    // translate to account for T0
    GPoint  trans{
       -(inv00*T0.x + inv01*T0.y),
       -(inv10*T0.x + inv11*T0.y)
    };

    GMatrix T_inv(col0, col1, trans);

    //  texture→device = triDev ∘ T_inv
    GMatrix texToDev = GMatrix::Concat(triDev, T_inv);

    // give to real shader
    return fReal->setContext(texToDev);
}



// delegate the actual sampling to the real shader
void ProxyShader::shadeRow(int x, int y, int count, GPixel row[]) {
fReal->shadeRow(x, y, count, row);
}


// --- ComposeShader Implementation ---
// --- ComposeShader Implementation ---
// --- ComposeShader Implementation ---

class ComposeShader : public GShader {
    public:
        ComposeShader(std::shared_ptr<GShader> a, std::shared_ptr<GShader> b)
            : fA(std::move(a)), fB(std::move(b)) {}
    
        bool isOpaque() override {
            return fA->isOpaque() && fB->isOpaque();
        }
    
        bool setContext(const GMatrix& ctm) override {
            return fA->setContext(ctm) && fB->setContext(ctm);
        }
    
        void shadeRow(int x, int y, int count, GPixel row[]) override {
            std::vector<GPixel> tmp(count);
            fA->shadeRow(x, y, count, row);
            fB->shadeRow(x, y, count, tmp.data());
    
            for (int i = 0; i < count; ++i) {
                GPixel a = row[i];
                GPixel b = tmp[i];
    
                int A = (GPixel_GetA(a) * GPixel_GetA(b) + 127) / 255;
                int R = (GPixel_GetR(a) * GPixel_GetR(b) + 127) / 255;
                int G = (GPixel_GetG(a) * GPixel_GetG(b) + 127) / 255;
                int B = (GPixel_GetB(a) * GPixel_GetB(b) + 127) / 255;
    
                row[i] = GPixel_PackARGB(A, R, G, B);
            }
        }
    
    private:
        std::shared_ptr<GShader> fA, fB;
    };
    

CombinedShader::CombinedShader(std::unique_ptr<GShader> c, std::unique_ptr<GShader> t)
: fColor(std::move(c)), fTex(std::move(t)) {}

bool CombinedShader::isOpaque() {
    return fColor->isOpaque() && fTex->isOpaque();
}

bool CombinedShader::setContext(const GMatrix& ctm) {
    // run both, even if one fails
    bool okColor = fColor->setContext(ctm);
    bool okTex   = fTex  ->setContext(ctm);
    return okColor && okTex;
}

void CombinedShader::shadeRow(int x, int y, int count, GPixel row[]) {
    // make sure temp buffers r right size
    fColorBuf.resize(count);
    fTexBuf  .resize(count);

    // fetch each child's premultiplied output
    fColor->shadeRow(x, y, count, fColorBuf.data());
    fTex  ->shadeRow(x, y, count, fTexBuf.data());
    // debug -- show first pixel of this scanline
    // if (count > 0) {
    //         auto col = fColorBuf[0], tex = fTexBuf[0];
    //         std::cerr << "[CombinedShader] y="<<y<<" x="<<x
    //                   <<"  color A,R,G,B=("
    //                   <<GPixel_GetA(col)<<","
    //                   <<GPixel_GetR(col)<<","
    //                   <<GPixel_GetG(col)<<","
    //                   <<GPixel_GetB(col)<<")"
    //                   <<"  tex A,R,G,B=("
    //                   <<GPixel_GetA(tex)<<","
    //                   <<GPixel_GetR(tex)<<","
    //                   <<GPixel_GetG(tex)<<","
    //                   <<GPixel_GetB(tex)<<")\n";
    //     }
    

    for (int i = 0; i < count; ++i) {
        // unpack 8-bit premultiplied channels
        int Ac = GPixel_GetA(fColorBuf[i]),
            Rc = GPixel_GetR(fColorBuf[i]),
            Gc = GPixel_GetG(fColorBuf[i]),
            Bc = GPixel_GetB(fColorBuf[i]);

        int At = GPixel_GetA(fTexBuf[i]),
            Rt = GPixel_GetR(fTexBuf[i]),
            Gt = GPixel_GetG(fTexBuf[i]),
            Bt = GPixel_GetB(fTexBuf[i]);

        int Aout = (Ac * At + 127) / 255;
        int Rout = (Rc * Rt + 127) / 255;
        int Gout = (Gc * Gt + 127) / 255;
        int Bout = (Bc * Bt + 127) / 255;

        //pack
        row[i] = GPixel_PackARGB(Aout, Rout, Gout, Bout);
    }
}




// --- makeTriangleShader factory function ---
// --- makeTriangleShader factory function ---
// --- makeTriangleShader factory function ---

std::unique_ptr<GShader> makeTriangleShader(
    const GPoint& P0, const GPoint& P1, const GPoint& P2,
    const GColor& C0, const GColor& C1, const GColor& C2,
    const GPoint& T0, const GPoint& T1, const GPoint& T2,
    GShader*   baseShader)
{
    const bool hasColor = (C0.a != 0) || (C1.a != 0) || (C2.a != 0);

    // only treat it as a “texture” case if:
    //  1) we were handed a shader, AND
    //  2) the three UVs are *not* all the same point
    const bool uvDegenerate = (T0.x == T1.x && T0.y == T1.y)
                           && (T0.x == T2.x && T0.y == T2.y);
    const bool hasTex      = (baseShader != nullptr) && !uvDegenerate;

    // 1) color only
    if (hasColor && !hasTex) {
        return std::make_unique<TriColorShader>(P0, P1, P2, C0, C1, C2);
    }

    // 2) texture only
    if (!hasColor && hasTex) {
        return std::make_unique<ProxyShader>(P0, P1, P2, T0, T1, T2, baseShader);
    }

    // 3) both color AND texture
    if (hasColor && hasTex) {
        // debug
        // std::cerr << "[makeTriangleShader] BOTH branch:\n"
        //           << " P0=(" << P0.x << "," << P0.y << ") C0=("
        //           << C0.r << "," << C0.g << "," << C0.b << ") T0=("
        //           << T0.x << "," << T0.y << ")\n"
        //           << " P1=(" << P1.x << "," << P1.y << ") C1=("
        //           << C1.r << "," << C1.g << "," << C1.b << ") T1=("
        //           << T1.x << "," << T1.y << ")\n"
        //           << " P2=(" << P2.x << "," << P2.y << ") C2=("
        //           << C2.r << "," << C2.g << "," << C2.b << ") T2=("
        //           << T2.x << "," << T2.y << ")\n";
        auto colorUP = std::unique_ptr<GShader>(
            new TriColorShader(P0, P1, P2, C0, C1, C2)
        );
        auto texUP = std::unique_ptr<GShader>(
            new ProxyShader(P0, P1, P2, T0, T1, T2, baseShader)
        );
        return std::make_unique<CombinedShader>(std::move(colorUP),
                                                std::move(texUP));
    }

    // 4) neither → nothing
    return nullptr;
}
