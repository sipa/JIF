#ifndef _BITSH_H_
#define _BITSH_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"


class TransformBS : public Transform {
public:
    bool virtual init(const ColorRanges *srcRanges) {
        if (srcRanges->numPlanes() < 3) return false;
        if (srcRanges->min(0) < 0 || srcRanges->min(1) < 0 || srcRanges->min(2) < 0) return false;
        return true;
    }

    void data(Image& image) const {
        printf("TransformBitShuffle::data\n");
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int R=image(0,r,c), G=image(1,r,c), B=image(2,r,c);
                // RRRRRRRR GGGGGGGG BBBBBBBB
                // GGRBGRBR GGRBGRBR GGRBBRBB
                image(0,r,c) = (G&128)+(G&64)+((R&128)>>2)+((B&128)>>3)+((G&32)>>2)+((R&64)>>4)+((B&64)>>5)+((R&32)>>5);
                R <<= 3;
                G <<= 3;
                B <<= 2;
                image(1,r,c) = (G&128)+(G&64)+((R&128)>>2)+((B&128)>>3)+((G&32)>>2)+((R&64)>>4)+((B&64)>>5)+((R&32)>>5);
                R <<= 3;
                G <<= 3;
                B <<= 2;
                image(2,r,c) = (G&128)+(G&64)+((R&128)>>2)+((B&128)>>3)+((B&64)>>3)+((R&64)>>4)+((B&32)>>4)+((B&16)>>4);
            }
        }
    }

    void invData(Image& image) const {
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                // GGRBGRBR GGRBGRBR GGRBBRBB
                // RRRRRRRR GGGGGGGG BBBBBBBB
                int X=image(0,r,c), Y=image(1,r,c), Z=image(2,r,c);
                image(0,r,c) = ((X&32)<<2)+((X&4)<<4)+((X&1)<<5)+((Y&32)>>1)+((Y&4)<<1)+((Y&1)<<2)+((Z&32)>>4)+((Z&4)>>2);
                image(1,r,c) = (X&128)+(X&64)+((X&8)<<2)+((Y&128)>>3)+((Y&64)>>3)+((Y&8)>>1)+((Z&128)>>6)+((Z&64)>>6);
                image(2,r,c) = ((X&16)<<3)+((X&2)<<5)+((Y&16)<<1)+((Y&2)<<3)+((Z&16)>>1)+((Z&8)>>1)+((Z&2))+((Z&1));
            }
        }
    }
};


#endif
