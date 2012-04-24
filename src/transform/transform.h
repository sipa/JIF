#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "../image/image.h"
#include "../image/color_range.h"
#include "../maniac/rac.h"

typedef RacInput40 RacIn;
typedef RacOutput40 RacOut;

class Transform {
protected:

public:
    virtual ~Transform() {};

    bool virtual init(const ColorRanges *srcRanges) { return true; }
    bool virtual process(const ColorRanges *srcRanges, const Image &image) { return true; };
    void virtual load(const ColorRanges *srcRanges, RacIn &rac) {};
    void virtual save(const ColorRanges *srcRanges, RacOut &rac) const {};
    const ColorRanges virtual *meta(Image& image, const ColorRanges *srcRanges) { return new DupColorRanges(srcRanges); }
    void virtual data(Image& image) const {}
    void virtual invData(Image& image) const {}

    // On save: init, process, save, meta, data, <processing>
    // On load: init, load, meta, <processing>, invData (reverse order)
};

#endif
