
/*
 *  Copyright 2025 Sarah Ouda
 */


#ifndef MYGPATH_DEFINED
#define MYGPATH_DEFINED

#include "include/GPath.h"  
#include "include/GMatrix.h"
#include <memory>
#include <vector>
#include <optional>


class MyGPath : public std::enable_shared_from_this<MyGPath> {
public:
    MyGPath(std::vector<GPoint> pts, std::vector<GPathVerb> vbs)
        : fPts(std::move(pts)), fVbs(std::move(vbs)) {}

    // returns bounds of path
    GRect bounds() const;

    // returns transformed path
    std::shared_ptr<MyGPath> transform(const GMatrix&) const;

    // helpers for subdividing curves
    static void ChopQuadAt(const GPoint src[3], GPoint dst[5], float t);
    static void ChopCubicAt(const GPoint src[4], GPoint dst[7], float t);

    // nested edger class for iterating over path edges
    class Edger {
    public:
        Edger(const MyGPath& path);
        std::optional<GPathVerb> next(GPoint pts[]);
    private:
        const GPoint*    fPrevMove;
        const GPoint*    fCurrPt;
        const GPathVerb* fCurrVb;
        const GPathVerb* fStopVb;
        int fPrevVerb;  // using -1 to indicate no previous verb
    };

    std::vector<GPoint> fPts;
    std::vector<GPathVerb> fVbs;
};

#endif
