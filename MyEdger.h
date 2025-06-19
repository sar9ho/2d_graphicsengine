/*
 *  Copyright 2025 Sarah Ouda
 */


#ifndef MYEDGER_DEFINED
#define MYEDGER_DEFINED

#include "include/GPath.h"
#include <vector>
#include <optional>

/**
 * uses prof. Reed's GPath::Edger to iterate over the path,
 * then post-processes any quadratic or cubic segments to flatten them into
 * line segments (with a 1/4-pixel tolerance)
 */

class MyEdger {
public:
    // construct the MyEdger from given GPath
    MyEdger(const GPath& path);

    // return next edge (always as a kLine segment) in pts (2 pnts)
    std::optional<GPathVerb> next(GPoint pts[]);
    
private:
    // instance of prof's Edger
    GPath::Edger fEdger;
    
    // buffer for flattened curve points, if curve was encountered
    std::vector<GPoint> fFlattenedPts;
    int fCurrentFlatIndex;
    GPoint fContourStart;
    GPoint fLastPt;
    bool fHasPendingClose = false;
};

#endif
