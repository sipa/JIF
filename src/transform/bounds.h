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
    void minmax(const int p, const prevPlanes &pp, ColorVal &min, ColorVal &max) const {
        if (p==0) { min=bounds[p].first; max=bounds[p].second; return; } // optimization for special case
        ranges->minmax(p, pp, min, max);
        if (min < bounds[p].first) min=bounds[p].first;
        if (max > bounds[p].second) max=bounds[p].second;
    }
/*
    ColorVal min(int p, int r, int c) const { return std::max(ranges->min(p,r,c), bounds[p].first); }
    ColorVal max(int p, int r, int c) const { return std::min(ranges->max(p,r,c), bounds[p].second); }
    ColorVal min(int p, int z, int r, int c) const { if (p==0) return bounds[p].first; else return std::max(ranges->min(p,z,r,c), bounds[p].first); }
    ColorVal max(int p, int z, int r, int c) const { if (p==0) return bounds[p].second; else return std::min(ranges->max(p,z,r,c), bounds[p].second); }
*/
//    ColorVal min(int p, int z, int r, int c) const { return min(p); }
//    ColorVal max(int p, int z, int r, int c) const { return max(p); }
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
        SimpleSymbolCoder<StaticBitChance, RacIn, 24> coder(rac);
        bounds.clear();
        for (int p=0; p<srcRanges->numPlanes(); p++) {
            ColorVal min = coder.read_int(0, srcRanges->max(p) - srcRanges->min(p)) + srcRanges->min(p);
            ColorVal max = coder.read_int(0, srcRanges->max(p) - min) + min;
            bounds.push_back(std::make_pair(min,max));
            fprintf(stdout,"plane[%i] : %i..%i \t",p,min,max);
        }
        fprintf(stdout,"\n");
    }

    void save(const ColorRanges *srcRanges, RacOut &rac) const {
        SimpleSymbolCoder<StaticBitChance, RacOut, 24> coder(rac);
        for (int p=0; p<srcRanges->numPlanes(); p++) {
            ColorVal min = bounds[p].first;
            ColorVal max = bounds[p].second;
            coder.write_int(0, srcRanges->max(p) - srcRanges->min(p), min - srcRanges->min(p));
            coder.write_int(0, srcRanges->max(p) - min, max - min);
            fprintf(stdout,"plane[%i] : %i..%i \t",p,min,max);
        }
        fprintf(stdout,"\n");
    }

    bool process(const ColorRanges *srcRanges, const Image &image) {
        bounds.clear();
        for (int p=0; p<srcRanges->numPlanes(); p++) {
            ColorVal min = srcRanges->max(p);
            ColorVal max = srcRanges->min(p);
            for (int r=0; r<image.rows(); r++) {
                for (int c=0; c<image.cols(); c++) {
#ifdef SMOOTHZOOM
                    ColorVal v = (image(p,r,c) & PARITYMASK);
#else
                    ColorVal v = image(p,r,c);
#endif
                    if (v < min) min = v;
                    if (v > max) max = v;
                    assert(v <= srcRanges->max(p));
                    assert(v >= srcRanges->min(p));
                }
            }
            bounds.push_back(std::make_pair(min,max));
        }
        return true;
    }
};

#endif
