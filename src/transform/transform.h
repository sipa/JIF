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

    void virtual load(const ColorRanges *srcRanges, RacIn &rac) {};
    void virtual process(const ColorRanges *srcRanges, const Image &image) {};
    void virtual save(const ColorRanges *srcRanges, RacOut &rac) const {};

public:
    virtual ~Transform() {};

    void virtual data(Image& image) const {}
    void virtual invData(Image& image) const {}

    bool virtual initFromImage(Image &image, RacOut& rac, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        process(srcRanges, image);
        save(srcRanges, rac);
        return meta(image, srcRanges, dstRanges);
    }

    bool virtual initFromRac(Image &image, RacIn& rac, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        load(srcRanges, rac);
        return meta(image, srcRanges, dstRanges);
    }
};

#endif
