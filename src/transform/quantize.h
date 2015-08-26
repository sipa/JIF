#ifndef _QUANT_H_
#define _QUANT_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"


class ColorRangesQuant : public ColorRanges
{
protected:
    const Image *image;
    const ColorRanges *ranges;
    int orig_numPlanes;
public:
    ColorRangesQuant(Image &imageIn, const ColorRanges *rangesIn) : image(&imageIn), ranges(rangesIn), orig_numPlanes(rangesIn->numPlanes()) {}
    bool isStatic() const { return false; }
    int numPlanes() const { return orig_numPlanes*2; }

    ColorVal min(int p) const { if (p<orig_numPlanes) return ranges->min(p)/QUANTIZATION;
                                else return 0; }
    ColorVal max(int p) const { if (p<orig_numPlanes) return ranges->max(p)/QUANTIZATION;
                                else return QUANTIZATION-1; }
};


class TransformQuantize : public Transform {
protected:
    int np;
public:
    bool virtual init(const ColorRanges *srcRanges) {
        np = srcRanges->numPlanes();
        return true;
    }

    const ColorRanges *meta(Image& image, const ColorRanges *srcRanges) {
        for (int p=0; p<np; p++) {
          image.add_plane(0,QUANTIZATION-1);
        }
        return new ColorRangesQuant(image, srcRanges);
    }

    void data(Image& image) const {
        printf("TransformQuantize::data: q=%i\n", QUANTIZATION);
        for (int p=0; p<np; p++) {
          for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                ColorVal pixel=image(p,r,c);
                image(p,r,c) = pixel >> QUANTIZATIONSHIFT;
                image(np+p,r,c) = pixel & QUANTIZATIONMASK;
            }
          }
        }
    }

    void invData(Image& image) const {
        for (int p=0; p<np; p++) {
          for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                ColorVal pixelA=image(p,r,c);
                ColorVal pixelB=image(np+p,r,c);
                image(p,r,c) = (pixelA << QUANTIZATIONSHIFT)+pixelB;
            }
          }
        }
        image.drop_planes(np);
    }
};


#endif
