#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "../image/image.h"
#include "../image/color_range.h"
#include "../maniac/rac.h"

typedef RacInput40 RacIn;
typedef RacOutput40 RacOut;

class Transform {
protected:
    bool virtual meta(Image& image, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        dstRanges = new DupColorRanges(srcRanges);
        return true;
    }

public:
    virtual ~Transform() {};

    void virtual data(Image& image) {}
    void virtual invData(Image& image) {}

    bool virtual initFromImage(Image &image, RacOut& rac, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        return meta(image, srcRanges, dstRanges);
    }

    bool virtual initFromRac(Image &image, RacIn& rac, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        return meta(image, srcRanges, dstRanges);
    }
};

#endif
