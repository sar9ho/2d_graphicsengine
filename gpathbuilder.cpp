/*
 *  Copyright 2025 Sarah Ouda
 */

 #include "include/GPathBuilder.h"
 #include "include/GMatrix.h"
 #include "include/GPath.h"
 #include <cassert>

 void GPathBuilder::addRect(const GRect& rect, GPathDirection dir) {
     // new contour must start at the top-left corner
     // For kCW, we traverse: top-left -> top-right -> bottom-right -> bottom-left
     // For kCCW, traverse: top-left -> bottom-left -> bottom-right -> top-right
     moveTo({ rect.left, rect.top });
     
     if (dir == GPathDirection::kCW) {
         lineTo({ rect.right, rect.top });
         lineTo({ rect.right, rect.bottom });
         lineTo({ rect.left, rect.bottom });
     } else {  // kCCW
         lineTo({ rect.left, rect.bottom });
         lineTo({ rect.right, rect.bottom });
         lineTo({ rect.right, rect.top });
     }
     // do not add an "close" command here
     // The drawing code (via the Edger) will auto-close the contour if needed
 }
 
 void GPathBuilder::addPolygon(const GPoint pts[], int count) {
     // Ensure there is at least one point
     if (count <= 0) {
         return;
     }
     // Start new contour with the first point
     moveTo(pts[0]);
     // Connect each subsequent point
     for (int i = 1; i < count; ++i) {
         lineTo(pts[i]);
     }
     // Again, auto-closing is handled later
 }
 
 void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection dir) {
    const float kappa = 0.5522847498f;

    float cx = center.x;
    float cy = center.y;
    float r = radius;
    float c = kappa * r;

    if (dir == GPathDirection::kCW) {
        moveTo({cx, cy - r});  // Top
        cubicTo({cx + c, cy - r}, {cx + r, cy - c}, {cx + r, cy});  // Top-right
        cubicTo({cx + r, cy + c}, {cx + c, cy + r}, {cx, cy + r});  // Bottom-right
        cubicTo({cx - c, cy + r}, {cx - r, cy + c}, {cx - r, cy});  // Bottom-left
        cubicTo({cx - r, cy - c}, {cx - c, cy - r}, {cx, cy - r});  // Top-left
    } else {
        moveTo({cx, cy - r});  // Top
        cubicTo({cx - c, cy - r}, {cx - r, cy - c}, {cx - r, cy});  // Top-left
        cubicTo({cx - r, cy + c}, {cx - c, cy + r}, {cx, cy + r});  // Bottom-left
        cubicTo({cx + c, cy + r}, {cx + r, cy + c}, {cx + r, cy});  // Bottom-right
        cubicTo({cx + r, cy - c}, {cx + c, cy - r}, {cx, cy - r});  // Top-right
    }
}


