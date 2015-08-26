#ifndef _SG_H_
#define _SG_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"


class TransformSG : public Transform {
public:
    bool virtual init(const ColorRanges *srcRanges) {
        if (srcRanges->numPlanes() < 3) return false;
        if (srcRanges->min(0) < 0 || srcRanges->min(1) < 0 || srcRanges->min(2) < 0) return false;
        return true;
    }

    void data(Image& image) const {
        printf("TransformSubtractGreen::data\n");
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int R=image(0,r,c), G=image(1,r,c), B=image(2,r,c);
                image(0,r,c) = (R-G) &0xFF;
                image(1,r,c) = G;
                image(2,r,c) = (B-G) &0xFF;
            }
        }
    }

    void invData(Image& image) const {
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int R=image(0,r,c), G=image(1,r,c), B=image(2,r,c);
                image(0,r,c) = (R+G) &0xFF;
                image(1,r,c) = G;
                image(2,r,c) = (B+G) &0xFF;
            }
        }
    }
};


#endif
