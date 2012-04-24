#ifndef _BOUNDS_H_
#define _BOUNDS_H_ 1

#include <vector>

#include "transform.h"
#include "../maniac/symbol.h"

class ColorRangesBounds : public ColorRanges
{
protected:
    std::vector<std::pair<ColorVal, ColorVal> > bounds;
    const ColorRanges *ranges;
public:
    ColorRangesBounds(const std::vector<std::pair<ColorVal, ColorVal> > &boundsIn, const ColorRanges *rangesIn) : bounds(boundsIn), ranges(rangesIn) {}
    bool isStatic() const { return false; }
    int numPlanes() const { return bounds.size(); }
    ColorVal min(int p) const { return std::max(ranges->min(p), bounds[p].first); }
    ColorVal max(int p) const { return std::min(ranges->max(p), bounds[p].second); }
    ColorVal min(int p, int r, int c) const { return std::max(ranges->min(p,r,c), bounds[p].first); }
    ColorVal max(int p, int r, int c) const { return std::min(ranges->max(p,r,c), bounds[p].second); }
};


class TransformBounds : public Transform {
protected:
    std::vector<std::pair<ColorVal, ColorVal> > bounds;

    const ColorRanges *meta(Image& image, const ColorRanges *srcRanges) {
        if (srcRanges->isStatic()) {
            return new StaticColorRanges(bounds);
        } else {
            return new ColorRangesBounds(bounds, srcRanges);
        }
    }

    void load(const ColorRanges *srcRanges, RacIn &rac) {
        SimpleSymbolCoder<StaticBitChance, RacIn> coder(rac, 24);
        bounds.clear();
        for (int p=0; p<srcRanges->numPlanes(); p++) {
            ColorVal min = coder.read_int(0, srcRanges->max(p) - srcRanges->min(p)) + srcRanges->min(p);
            ColorVal max = coder.read_int(0, srcRanges->max(p) - min) + min;
            bounds.push_back(std::make_pair(min,max));
            fprintf(stderr,"plane[%i] : %i..%i\n",p,min,max);
        }
    }

    void save(const ColorRanges *srcRanges, RacOut &rac) const {
        SimpleSymbolCoder<StaticBitChance, RacOut> coder(rac, 24);
        for (int p=0; p<srcRanges->numPlanes(); p++) {
            ColorVal min = bounds[p].first;
            ColorVal max = bounds[p].second;
            coder.write_int(0, srcRanges->max(p) - srcRanges->min(p), min - srcRanges->min(p));
            coder.write_int(0, srcRanges->max(p) - min, max - min);
            fprintf(stderr,"plane[%i] : %i..%i\n",p,min,max);
        }
    }

    bool process(const ColorRanges *srcRanges, const Image &image) {
        bounds.clear();
        for (int p=0; p<srcRanges->numPlanes(); p++) {
            ColorVal min = srcRanges->max(p);
            ColorVal max = srcRanges->min(p);
            for (int r=0; r<image.rows(); r++) {
                for (int c=0; c<image.cols(); c++) {
                    ColorVal v = image(p,r,c);
                    if (v < min) min = v;
                    if (v > max) max = v;
                }
            }
            bounds.push_back(std::make_pair(min,max));
        }
        return true;
    }
};

#endif
