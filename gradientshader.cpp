/*
 *  Copyright 2025 Sarah Ouda
 */

#include "include/GShader.h"
#include "include/GMatrix.h"
#include "include/GPoint.h"
#include "include/GColor.h"
#include "include/GPixel.h"
#include "GPixelUtils.h"
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>
#include <optional>


// implements gradient shader that interpolates between colors along a line defined by two endpoints
class LinearGradientShader : public GShader {
public:
    // Default to kClamp if no tile mode is passed
    LinearGradientShader(GPoint p0, GPoint p1, const GColor colors[], int count)
        : LinearGradientShader(p0, p1, colors, count, GTileMode::kClamp) {}

    // Full constructor
    LinearGradientShader(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode)
        : fP0(p0), fP1(p1), fColors(colors, colors + count), fCount(count), fTileMode(mode) {


    float dx = fP1.x - fP0.x;
    float dy = fP1.y - fP0.y;
    float L = std::sqrt(dx*dx + dy*dy);
    if (L == 0) {
        fLocalMatrix = GMatrix();
    } else {
        GMatrix T = GMatrix::Translate(-fP0.x, -fP0.y);
        float theta = std::atan2(dy, dx);
        GMatrix R = GMatrix::Rotate(-theta);
        GMatrix S = GMatrix::Scale(1.0f / L, 1.0f / L);
        fLocalMatrix = GMatrix::Concat(S, GMatrix::Concat(R, T));
    }

    fValid = false;
    }


    virtual bool isOpaque() override {
        // The gradient is opaque if every color has alpha == 1
        for (const auto& c : fColors) {
            if (c.a < 1.0f) {
                return false;
            }
        }
        return true;
    }

    virtual bool setContext(const GMatrix& ctm) override {
        // Combine the current CTM with our local matrix----
        //  maps points from our gradient space (unit space) into device space
        //  then invert the combined matrix to map device coordinates into gradient space
        std::optional<GMatrix> ctmInv = ctm.invert();
        if (!ctmInv.has_value()) {
            return false;
        }
        fCombined = GMatrix::Concat(fLocalMatrix, ctmInv.value());
        fValid = true;

        return true;
    }

    virtual void shadeRow(int x, int y, int count, GPixel row[]) override {
        if (!fValid) {
            for (int i = 0; i < count; ++i) {
                row[i] = 0;
            }
            return;
        }
        // Map the leftmost pixel's center from device space into gradient space
        GPoint devicePt = { float(x) + 0.5f, float(y) + 0.5f };
        GPoint localPt;
        fCombined.mapPoints(&localPt, &devicePt, 1);

        // For a properly aligned gradient, we always use the x-coordinate
        float t = localPt.x;
        
        // Compute dt by mapping the difference from one pixel over
        GPoint devicePtNext = { float(x + 1) + 0.5f, float(y) + 0.5f };
        GPoint localPtNext;
        fCombined.mapPoints(&localPtNext, &devicePtNext, 1);
        float dt = localPtNext.x - localPt.x;
     
        for (int i = 0; i < count; ++i) {
            float finalT;
            switch(fTileMode) {
                case GTileMode::kClamp:
                    finalT = std::max(0.0f, std::min(1.0f, t));
                    break;
                case GTileMode::kRepeat:
                    finalT = t - std::floor(t);
                    break;
                case GTileMode::kMirror: {
                    float mod = std::fmod(t, 2.0f);
                    if (mod < 0) { 
                        mod += 2.0f; 
                    }
                    finalT = (mod > 1.0f) ? 2.0f - mod : mod;
                    break;
                }
            }
        
            GColor color;
            if (fCount == 1) {
                color = fColors[0];
            } else {
                float gap = 1.0f / (fCount - 1);
                int idx = std::clamp(int(finalT * (fCount - 1)), 0, fCount - 2);
                float tLocal = (finalT - idx * gap) / gap;
                color = lerpColor(fColors[idx], fColors[idx + 1], tLocal);
            }
            row[i] = MultandScale(color);
            t += dt;
        }
        
    }
    

private:
    GPoint fP0, fP1;
    std::vector<GColor> fColors;
    int fCount;
    GMatrix fLocalMatrix;   // Maps world coordinates to unit gradient space (p0 -> (0,0), p1 -> (1,0))
    GMatrix fCombined;      // Inverse of (CTM concat fLocalMatrix), used to map device points to gradient space
    bool fValid;
    GTileMode fTileMode;

    // Helper for linear interpolation between two GColors
    GColor lerpColor(const GColor& c0, const GColor& c1, float t) {
        return {
            c0.r * (1 - t) + c1.r * t,
            c0.g * (1 - t) + c1.g * t,
            c0.b * (1 - t) + c1.b * t,
            c0.a * (1 - t) + c1.a * t,
        };
    }
};

// Factory function --- declared in GShader.h
std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count) {
    if (count < 1) {
        return nullptr;
    }
    return std::make_shared<LinearGradientShader>(p0, p1, colors, count);
}


std::shared_ptr<GShader> GCreateLinearGradient(GPoint p0, GPoint p1, const GColor colors[], int count, GTileMode mode) {
    if (count < 1) return nullptr;
    return std::make_shared<LinearGradientShader>(p0, p1, colors, count, mode);
}

