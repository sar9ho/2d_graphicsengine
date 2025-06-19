/*
 *  Copyright 2025 Sarah Ouda
 */


#include "GPixelUtils.h"
#include <cmath>      
#include <algorithm>  

GPixel MultandScale(const GColor& color) {
    float r = color.r;
    float g = color.g;
    float b = color.b;
    float a = color.a;

    int r_premult = GRoundToInt(r * a * 255);
    int g_premult = GRoundToInt(g * a * 255);
    int b_premult = GRoundToInt(b * a * 255);
    int alpha     = GRoundToInt(a * 255);

    return GPixel_PackARGB(alpha, r_premult, g_premult, b_premult);
}
