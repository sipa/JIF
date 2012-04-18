#ifndef _YIQ_H_
#define _YIQ_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"

class ColorRangesYIQ : public ColorRanges
{
protected:
    const Image *image;
public:
    ColorRangesYIQ(Image &imageIn) : image(&imageIn) {}
    bool isStatic() const { return false; }
    ColorVal min(int p) const { return 0; }
    ColorVal max(int p) const { switch(p) { case 0: return 255; case 1: return 510; case 2: return 510; }; assert(false); return 0; }
    ColorVal min(int p, int r, int c) const;
    ColorVal max(int p, int r, int c) const;
};


class TransformYIQ : public Transform {
protected:
    bool data(Image& image) {
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int R=image(0,r,c), G=image(1,r,c), B=image(2,r,c);
                int Y = ((R + B) / 2 + G) / 2;
                int I = R - B + 255;
                int Q = (R + B) / 2 - G + 255;
                image(0,r,c) = Y;
                image(1,r,c) = I;
                image(2,r,c) = Q;
            }
        }
        return true;
    }

public:
//    bool full(Image& image, const ColorRanges *srcRanges, const ColorRanges& *dstRanges) {
//        return (meta(srcRanges, image, dstRanges) && data(image));
//    }

    bool meta(Image& image, const ColorRanges *srcRanges, const ColorRanges *&dstRanges) {
        if (image.numPlanes() != 3) return false;
        if (srcRanges->min(0) < 0 || srcRanges->max(0) > 255) return false;
        if (srcRanges->min(1) < 0 || srcRanges->max(1) > 255) return false;
        if (srcRanges->min(2) < 0 || srcRanges->max(2) > 255) return false;

        dstRanges = new ColorRangesYIQ(image);
        return true;
    }

    bool invData(Image& image) {
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int Y=image(0,r,c), I=image(1,r,c), Q=image(2,r,c);
                int R = Y + (Q + 2) / 2 + (I + 2) / 2 - 256;
                int G = Y - (Q + 1) / 2 + 128;
                int B = Y + (Q + 2) / 2 - (I + 1) / 2;
                image(0,r,c) = R;
                image(1,r,c) = G;
                image(2,r,c) = B;
            }
        }
        return true;
    }
};


#endif
