#ifndef _PP_H_
#define _PP_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"


class ColorRangesPP : public ColorRanges
{
protected:
    std::vector<int> permutation;
    const ColorRanges *ranges;
public:
    ColorRangesPP(const std::vector<int> permIn, const ColorRanges *rangesIn) : permutation(permIn), ranges(rangesIn) {}
    bool isStatic() const { return false; }
    int numPlanes() const { return permutation.size(); }
    ColorVal min(int p) const { return ranges->min(permutation[p]); }
    ColorVal max(int p) const { return ranges->max(permutation[p]);  }
    ColorVal min(int p, int r, int c) const { return ranges->min(permutation[p]); }
    ColorVal max(int p, int r, int c) const { return ranges->max(permutation[p]); }
    ColorVal min(int p, int z, int r, int c) const { return ranges->min(permutation[p]); }
    ColorVal max(int p, int z, int r, int c) const { return ranges->max(permutation[p]); }
};


class TransformPP : public Transform {
protected:
    std::vector<int> permutation;
    std::vector<int> inv_permutation;

public:
    bool virtual init(const ColorRanges *srcRanges) {
        if (srcRanges->numPlanes() != 3) return false;
        permutation.push_back(0);  // new 0 = old A
        permutation.push_back(1);  // new 1 = old B
        permutation.push_back(2);  // new 2 = old C

        inv_permutation.push_back(0);  // old 0 = new A
        inv_permutation.push_back(1);  // old 1 = new B
        inv_permutation.push_back(2);  // old 2 = new C
        return true;
    }

    const ColorRanges *meta(Image& image, const ColorRanges *srcRanges) {
        return new ColorRangesPP(permutation, srcRanges);
    }

    void data(Image& image) const {
        printf("Transform Permute Planes\n");
        image.permute_planes(permutation);
    }

    void invData(Image& image) const {
        printf("Transform Permute Planes\n");
        image.permute_planes(inv_permutation);
    }
};


#endif
