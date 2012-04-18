#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "../image/image.h"
#include "../image/color_range.h"

class Transform {
protected:
    bool virtual data(Image& image) {
        return true;
    }

public:
    virtual ~Transform() {};

    bool virtual full(Image& image, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        return (meta(image, srcRanges, dstRanges) && data(image));
    }

    bool virtual meta(Image& image, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        dstRanges = new DupColorRanges(srcRanges);
        return true;
    }

    bool virtual invData(Image& image) {
        return true;
    }
};

#endif
