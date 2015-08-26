#ifndef _YIQ_H_
#define _YIQ_H_ 1

#include "../image/image.h"
#include "../image/color_range.h"
#include "transform.h"

#define clip(x,l,u)   if (x < l) x=l; if (x > u) x=u

ColorVal static inline get_min_y(int par) {
    return 0;
}

ColorVal static inline get_max_y(int par) {
    return par*4-1;
}

ColorVal static inline get_min_i(int par, ColorVal y) {
    assert(y >= get_min_y(par));
    assert(y <= get_max_y(par));

    if (y<par-1) {
      return 4*par-4-4*y;
    } else if (y>=3*par) {
      return 3+4*(y-3*par);
    } else {
      return 0;
    }
}

ColorVal static inline get_max_i(int par, ColorVal y) {
    assert(y >= get_min_y(par));
    assert(y <= get_max_y(par));

    if (y<par-1) {
      return 4*par+2+4*y;
    } else if (y>=3*par) {
      return 8*par-5-4*(y-3*par);
    } else {
      return 8*par-2;
    }
}

ColorVal static inline get_min_q(int par, ColorVal y, ColorVal i) {
    assert(y >= get_min_y(par));
    assert(y <= get_max_y(par));
    if (i < get_min_i(par,y)) return 8*par; //invalid value
    if (i > get_max_i(par,y)) return 8*par; //invalid value
    assert(i >= get_min_i(par,y));
    assert(i <= get_max_i(par,y));

    if (y<par-1) {
      return 4*par-2-2*y+(abs(i-4*par+1)/2)*2;
    } else if (y>=3*par) {
      return 4*par-1-2*(4*par-1-y);
    } else {
      return std::max(1+(y-2*par)*2, 2*par-(y-par+1)*2+(abs(i-4*par+1)/2)*2);
    }
}

ColorVal static inline get_max_q(int par, ColorVal y, ColorVal i) {
    assert(y >= get_min_y(par));
    assert(y <= get_max_y(par));
    if (i < get_min_i(par,y)) return -1; //invalid value
    if (i > get_max_i(par,y)) return -1; //invalid value
    assert(i >= get_min_i(par,y));
    assert(i <= get_max_i(par,y));

    if (y<par-1) {
      return 4*par+2*y;
    } else if (y>=3*par) {
      return 4*par-1+2*(4*par-1-y)-((1+abs(i-4*par+1))/2)*2;
    } else {
      return std::min(6*par-2+(y-par+1)*2, 6*par-1+(3*par-1-y)*2-((1+abs(i-4*par+1))/2)*2);
    }
}


class ColorRangesYIQ : public ColorRanges
{
protected:
    const Image *image;
    const int par=64; // range: [0..4*par-1]
    const ColorRanges *ranges;
public:
//    ColorRangesYIQ(Image &imageIn, int parIn, const ColorRanges *rangesIn) : image(&imageIn), par(parIn), ranges(rangesIn) {}
    ColorRangesYIQ(Image &imageIn, int parIn, const ColorRanges *rangesIn) : image(&imageIn), ranges(rangesIn) { if (parIn != par) printf("OOPS: using YIQ transform on something other than rgb888 ?\n");}
    bool isStatic() const { return false; }
    int numPlanes() const { return ranges->numPlanes(); }

    ColorVal min(int p) const { if (p<3) return 0; else return ranges->min(p); }
    ColorVal max(int p) const { switch(p) {
                                        case 0: return 4*par-1;
                                        case 1: return 8*par-2;
                                        case 2: return 8*par-2;
                                        default: return ranges->max(p);
                                         };
                              }
    void minmax(const int p, const prevPlanes &pp, ColorVal &minv, ColorVal &maxv) const {
         if (p==1) { minv=get_min_i(par, pp[0]); maxv=get_max_i(par, pp[0]); return; }
         else if (p==2) { minv=get_min_q(par, pp[0], pp[1]); maxv=get_max_q(par, pp[0], pp[1]); return; }
         else if (p==0) { minv=0; maxv=get_max_y(par); return;}
         else ranges->minmax(p,pp,minv,maxv);
    }

/*
    ColorVal min(int p) const { return 0; }
    ColorVal max(int p) const { switch(p) { case 0: return 1024; case 1: return 2048; case 2: return 2048; }; assert(false); return 0; }
*/
/*
    ColorVal min(int p) const { return 0; }
    ColorVal max(int p) const { switch(p) { case 0: return 425; case 1: return 1020; case 2: return 1020; }; assert(false); return 0; }
*/
/*
#ifndef SMOOTHZOOM
#ifndef PERMUTEPLANES
    ColorVal min(int p, int r, int c) const;
    ColorVal max(int p, int r, int c) const;
    ColorVal min(int p, int z, int r, int c) const;
    ColorVal max(int p, int z, int r, int c) const;
#endif
#endif
*/
};


class TransformYIQ : public Transform {
protected:
    int par;

public:
    bool virtual init(const ColorRanges *srcRanges) {
        if (srcRanges->numPlanes() < 3) return false;
        if (srcRanges->min(0) < 0 || srcRanges->min(1) < 0 || srcRanges->min(2) < 0) return false;
        int max = std::max(std::max(srcRanges->max(0), srcRanges->max(1)), srcRanges->max(2));
        par = max/4+1;
        return true;
    }

    const ColorRanges *meta(Image& image, const ColorRanges *srcRanges) {
        return new ColorRangesYIQ(image, par, srcRanges);
    }

    void data(Image& image) const {
//        printf("TransformYIQ::data: par=%i\n", par);
        for (int r=0; r<image.rows(); r++) {
            for (int c=0; c<image.cols(); c++) {
                int R=image(0,r,c), G=image(1,r,c), B=image(2,r,c);

                int Y = ((R + B) / 2 + G) / 2;
                int I = R - B + par*4 - 1;
                int Q = (R + B) / 2 - G + par*4 - 1;

/*
                int Y = 0.5+4.0*(0.30*R + 0.59*G + 0.11*B);
                int I = 0.5+4.0*(1.00*R - 0.46*G - 0.54*B)+1024;
                int Q = 0.5+4.0*(0.40*R - 1.00*G + 0.60*B)+1024;
*/
/*
                int Y = (3*R  + 6*G + 1*B )/6;
                int I = (10*R - 5*G - 5*B )/5 + 510;
                int Q = (5*R  -10*G + 5*B )/5 + 510;
*/
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

/*
                float y = 0.25*Y;
                float i = 0.25*(I-1024);
                float q = 0.25*(Q-1024);
                int R = y + 0.568627*i + 0.328431*q +0.5;
                int G = y - 0.166667*i - 0.333333*q +0.5;
                int B = y - 0.656863*i + 0.892157*q +0.5;
*/
/*
                I -= 510;
                Q -= 510;
                Y *=6;
                I *=5;
                Q *=5;
                int R = 15*Y +  8*I +  5*Q +75;
                int G = 15*Y -  2*I -  5*Q +75;
                int B = 15*Y - 12*I + 15*Q +75;
                R /= 150;
                G /= 150;
                B /= 150;
*/
                // clipping only needed in case of lossy/partial decoding
                clip(R, 0, par*4-1);
                clip(G, 0, par*4-1);
                clip(B, 0, par*4-1);
                image(0,r,c) = R;
                image(1,r,c) = G;
                image(2,r,c) = B;
            }
        }
    }
};


#endif
