#ifndef _YIQ_H_
#define _YIQ_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"

class ColorRangesYIQ : public ColorRanges
{
protected:
    const Image *image;
    int par; // range: [0..4*par-1]
public:
    ColorRangesYIQ(Image &imageIn, int parIn) : image(&imageIn), par(parIn) {}
    bool isStatic() const { return false; }
    int numPlanes() const { return 3; }
    ColorVal min(int p) const { return 0; }
    ColorVal max(int p) const { switch(p) { case 0: return 4*par-1; case 1: return 8*par-2; case 2: return 8*par-2; }; assert(false); return 0; }
    ColorVal min(int p, int r, int c) const;
    ColorVal max(int p, int r, int c) const;
};


class TransformYIQ : public Transform {
protected:
    int par;

public:
    bool virtual init(const ColorRanges *srcRanges) {
        if (srcRanges->numPlanes() != 3) return false;
        if (srcRanges->min(0) < 0 || srcRanges->min(1) < 0 || srcRanges->min(2) < 0) return false;
        int max = std::max(std::max(srcRanges->max(0), srcRanges->max(1)), srcRanges->max(2));
        par = max/4+1;
        return true;
    }

    const ColorRanges *meta(Image& image, const ColorRanges *srcRanges) {
        return new ColorRangesYIQ(image, par);
    }

    void data(Image& image) const {
        printf("TransformYIQ::data: par=%i\n", par);
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int R=image(0,r,c), G=image(1,r,c), B=image(2,r,c);
                int Y = ((R + B) / 2 + G) / 2;
                int I = R - B + par*4 - 1;
                int Q = (R + B) / 2 - G + par*4 - 1;
                image(0,r,c) = Y;
                image(1,r,c) = I;
                image(2,r,c) = Q;
            }
        }
    }

    void invData(Image& image) const {
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int Y=image(0,r,c), I=image(1,r,c), Q=image(2,r,c);
                int R = Y + (Q + 2) / 2 + (I + 2) / 2 - 4*par;
                int G = Y - (Q + 1) / 2 + 2*par;
                int B = Y + (Q + 2) / 2 - (I + 1) / 2;
                image(0,r,c) = R;
                image(1,r,c) = G;
                image(2,r,c) = B;
            }
        }
    }
};


#endif
