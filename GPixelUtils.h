/*
 *  Copyright 2025 Sarah Ouda
 */

#ifndef GPIXELUTILS_H
#define GPIXELUTILS_H

#include "include/GColor.h"
#include "include/GPixel.h"

// Convert a GColor (unpremultiplied) to a premultiplied GPixel
GPixel MultandScale(const GColor& color);

#endif
