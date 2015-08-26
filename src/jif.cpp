#include <string>

#include "maniac/rac.h"
#include "maniac/compound.h"
#include "maniac/util.h"

#include "image/image.h"
#include "image/color_range.h"
#include "transform/factory.h"

#include "jif_config.h"



FILE *f;  // the compressed file

std::vector<ColorVal> grey; // a pixel with values in the middle of the bounds

typedef SimpleBitChance                         JifBitChancePass1;

// faster:
//typedef SimpleBitChance                         JifBitChancePass2;
//typedef SimpleBitChance                         JifBitChanceParities;
//typedef SimpleBitChance                         JifBitChanceMeta;

// better compression:
//typedef MultiscaleBitChance<6,SimpleBitChance>  JifBitChancePass1;
typedef MultiscaleBitChance<6,SimpleBitChance>  JifBitChancePass2;
typedef MultiscaleBitChance<6,SimpleBitChance>  JifBitChanceParities;
typedef MultiscaleBitChance<6,SimpleBitChance>  JifBitChanceMeta;

typedef MultiscaleBitChance<6,SimpleBitChance>  JifBitChanceTree;

template<typename RAC> void static write_name(RAC& rac, std::string str)
{
    UniformSymbolCoder<RAC> coder(rac);
    coder.write_int(3, 8, str.size());
    for (unsigned int i=0; i<str.size(); i++) {
        char c = str[i];
        int n = ((c >= 'A' && c <= 'Z') ? c - 'A' :
                ((c >= 'a' && c <= 'z') ? c - 'a' :
                ((c >= '0' && c <= '9') ? c - '0' + 26 : 36)));
        coder.write_int(0, 36, n);
    }
}

template<typename RAC> std::string static read_name(RAC& rac)
{
    static char cs[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
    UniformSymbolCoder<RAC> coder(rac);
    int l = coder.read_int(3, 8);
    std::string str;
    for (int i=0; i<l; i++) {
        int n = coder.read_int(0, 36);
        str += cs[n];
    }
    return str;
}




/******************************************/
/*   FFV1 encoding/decoding               */
/******************************************/


const int NB_PROPERTIES_FFV1[] = {7,8,9,7};
const int NB_PROPERTIES_FFV1A[] = {8,9,10,7};
void static initPropRanges_ffv1(Ranges &propRanges, const ColorRanges &ranges, int p)
{
    propRanges.clear();
    int min = ranges.min(p);
    int max = ranges.max(p);
    int mind = min - max, maxd = max - min;

    if (p != 3) {
      for (int pp = 0; pp < p; pp++) {
        propRanges.push_back(std::make_pair(ranges.min(pp), ranges.max(pp)));  // pixels on previous planes
      }
      if (ranges.numPlanes()>3) propRanges.push_back(std::make_pair(ranges.min(3), ranges.max(3)));  // pixel on alpha plane
    }
    propRanges.push_back(std::make_pair(min,max));   // guess (median of 3)
    propRanges.push_back(std::make_pair(0,3));       // which predictor was it
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
}


template<typename I> void static swap(I& a, I& b)
{
    I c = a;
    a = b;
    b = c;
}

template<typename I> I static median3(I a, I b, I c)
{
    if (a<b) swap(a,b);
    if (b<c) swap(b,c);
    if (a<b) swap(a,b);
//    assert(b>=0);
    return b;
}
/*
ColorVal static predict_ffv1(const Image &image, int p, int r, int c)
{
    ColorVal left = image(p,r,c-1);
    ColorVal top = image(p,r-1,c);
    ColorVal gradient = left + top - image(p,r-1,c-1);
    return median3(left,top,gradient);
}
void static calcProps_origFFV1(Properties &properties, const Image &image, int p, int r, int c)
{
    for (int pp = 0; pp < p; pp++) {
        properties.push_back(image(pp,r,c));
    }
    properties.push_back(image(p,r,c-1)-image(p,r-1,c-1));  // left - topleft
    properties.push_back(image(p,r-1,c-1)-image(p,r-1,c));  // topleft - top
    properties.push_back(image(p,r-1,c)-image(p,r-1,c+1));  // top - topright
    properties.push_back(image(p,r-2,c)-image(p,r-1,c));    // toptop - top
    properties.push_back(image(p,r,c-2)-image(p,r,c-1));    // leftleft - left
}
*/
ColorVal predict_and_calcProps_ffv1(Properties &properties, const ColorRanges *ranges, const Image &image, const int p, const int r, const int c, ColorVal &min, ColorVal &max) {
    ColorVal guess;
    int which = 0;
//    properties.clear();
    int index=0;
    if (p != 3) {
      for (int pp = 0; pp < p; pp++) {
        properties[index++] = image(pp,r,c);
      }
      if (image.numPlanes()>3) properties[index++] = image(3,r,c);
    }
    ColorVal left = (c>0 ? image(p,r,c-1) : grey[p]);;
    ColorVal top = (r>0 ? image(p,r-1,c) : grey[p]);
    ColorVal topleft = (r>0 && c>0 ? image(p,r-1,c-1) : grey[p]);
    ColorVal gradientTL = left + top - topleft;
    guess = median3(gradientTL, left, top);
    ranges->snap(p,properties,min,max,guess);
    if (guess == gradientTL) which = 0;
    else if (guess == left) which = 1;
    else if (guess == top) which = 2;


    properties[index++] = guess;
    properties[index++] = which;

    if (c > 0 && r > 0) { properties[index++] = left - topleft; properties[index++] = topleft - top; }
                 else   { properties[index++] = 0; properties[index++] = 0;  }

    if (c+1 < image.cols() && r > 0) properties[index++] = top - image(p,r-1,c+1); // top - topright
                 else   properties[index++] = 0;
    if (r > 1) properties[index++] = image(p,r-2,c)-top;    // toptop - top
         else properties[index++] = 0;
    if (c > 1) properties[index++] = image(p,r,c-2)-left;    // leftleft - left
         else properties[index++] = 0;

    return guess;
}

template<typename Coder> void encode_ffv1_inner(std::vector<Coder*> &coders, const Image &image, const ColorRanges *ranges)
{
    ColorVal min,max;
    long fs = ftell(f);
    long pixels = image.cols()*image.rows();
    int nump = image.numPlanes();
    int beginp = (nump>3 ? 3 : 0); int i=0;
    for (int p = beginp; i++ < nump; p = (p+1)%nump) {
        Properties properties((nump>3?NB_PROPERTIES_FFV1A[p]:NB_PROPERTIES_FFV1[p]));
        if (ranges->min(p) < ranges->max(p))
          fprintf(stdout,"[%i] ENC_FFV1_STYLE ",p);
        else continue;
        fflush(stdout);
        for (int r = 0; r < image.rows(); r++) {
            for (int c = 0; c < image.cols(); c++) {
                if (nump>3 && p<3 && image(3,r,c) == 0) continue;
                ColorVal guess = predict_and_calcProps_ffv1(properties,ranges,image,p,r,c,min,max);
                ColorVal curr = image(p,r,c);
                coders[p]->write_int(properties, min - guess, max - guess, curr - guess);
            }
        }
        long nfs = ftell(f);
        if (nfs-fs > 0) fprintf(stdout,"\tfilesize : %li (+%li for %li pixels, %f bpp)\n", nfs, nfs-fs, pixels, 8.0*(nfs-fs)/pixels );
        fs = nfs;
    }
}

template<typename Rac, typename Coder> void encode_ffv1_pass(Rac &rac, const Image &image, const ColorRanges *ranges, std::vector<Tree> &forest, int repeats)
{
    std::vector<Coder*> coders;
    for (int p = 0; p < ranges->numPlanes(); p++) {
        Ranges propRanges;
        initPropRanges_ffv1(propRanges, *ranges, p);
        coders.push_back(new Coder(rac, propRanges, forest[p]));
    }

//    encode_ffv1_inner(coders, image, ranges);
    if (repeats>1) printf("Iterating %i times to find a better tree.\n",repeats);
    while(repeats-- > 0) {
     encode_ffv1_inner(coders, image, ranges);
    }

    for (int p = 0; p < image.numPlanes(); p++) {
        coders[p]->simplify();
    }

    for (int p = 0; p < image.numPlanes(); p++) {
#ifdef STATS
        indent(0); printf("Plane %i\n", p);
        coders[p]->info(0+1);
#endif
        delete coders[p];
    }
}




void encode_ffv1_interpol_zero_alpha(Image &image, const ColorRanges *ranges)
{

    ColorVal min,max;
    int nump = image.numPlanes();
    if (nump > 3)
    for (int p = 0; p < 3; p++) {
        Properties properties((nump>3?NB_PROPERTIES_FFV1A[p]:NB_PROPERTIES_FFV1[p]));
        if (ranges->min(p) < ranges->max(p))
          fprintf(stdout,"[%i] interpol_zero_alpha ",p);
        else continue;
        for (int r = 0; r < image.rows(); r++) {
            for (int c = 0; c < image.cols(); c++) {
                if (image(3,r,c) == 0) {
                    image(p,r,c) = predict_and_calcProps_ffv1(properties,ranges,image,p,r,c,min,max);
                }
            }
        }
    }
}


template<typename Coder> void decode_ffv1_inner(std::vector<Coder*> &coders, Image &image, const ColorRanges *ranges)
{

    ColorVal min,max;
    int nump = image.numPlanes();
    int beginp = (nump>3 ? 3 : 0); int i=0;
    for (int p = beginp; i++ < nump; p = (p+1)%nump) {
        Properties properties((nump>3?NB_PROPERTIES_FFV1A[p]:NB_PROPERTIES_FFV1[p]));
        if (ranges->min(p) < ranges->max(p))
          fprintf(stdout,"[%i] DEC_FFV1_STYLE ",p);
        else continue;
        for (int r = 0; r < image.rows(); r++) {
            for (int c = 0; c < image.cols(); c++) {
                ColorVal guess = predict_and_calcProps_ffv1(properties,ranges,image,p,r,c,min,max);
                if (nump>3 && p<3 && image(3,r,c) == 0) { image(p,r,c)=guess; continue;}
                ColorVal curr = coders[p]->read_int(properties, min - guess, max - guess) + guess;
                image(p,r,c) = curr;
            }
        }
    }
}

template<typename Rac, typename Coder> void decode_ffv1_pass(Rac &rac, Image &image, const ColorRanges *ranges, std::vector<Tree> &forest)
{
    std::vector<Coder*> coders;
    for (int p = 0; p < image.numPlanes(); p++) {
        Ranges propRanges;
        initPropRanges_ffv1(propRanges, *ranges, p);
        coders.push_back(new Coder(rac, propRanges, forest[p]));
    }

    decode_ffv1_inner(coders, image, ranges);

    for (int p = 0; p < image.numPlanes(); p++) {
        delete coders[p];
    }
}




/******************************************/
/*   JIF2 encoding/decoding               */
/******************************************/
const int NB_PROPERTIES[] = {8,7,8,8};
const int NB_PROPERTIESA[] = {9,8,9,8};
void static initPropRanges(Ranges &propRanges, const ColorRanges &ranges, int p)
{
    propRanges.clear();
    int min = ranges.min(p);
    int max = ranges.max(p);
    int mind = min - max, maxd = max - min;
//    printf("initpropranges(%i): %i..%i\n",p,min,max);

    if (p != 3) {       // alpha channel first
      for (int pp = 0; pp < p; pp++) {
        propRanges.push_back(std::make_pair(ranges.min(pp), ranges.max(pp)));  // pixels on previous planes
      }
      if (ranges.numPlanes()>3) propRanges.push_back(std::make_pair(ranges.min(3), ranges.max(3)));  // pixel on alpha plane
    }
/*
    for (int pp = 0; pp < p; pp++) {
        int minpp = ranges.min(pp);
        int maxpp = ranges.max(pp);
        int mindpp = minpp - maxpp, maxdpp = maxpp - minpp;
        propRanges.push_back(std::make_pair(mindpp,maxdpp));        // pixel diffs on previous planes
    }
*/
    propRanges.push_back(std::make_pair(mind,maxd)); // neighbor A - neighbor B   (top-bottom or left-right)
    propRanges.push_back(std::make_pair(min,max));   // guess (median of 3)
//    propRanges.push_back(std::make_pair(0,max-min));   // guess (median of 3)
    propRanges.push_back(std::make_pair(0,3));       // which predictor was it
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
//    propRanges.push_back(std::make_pair(mind,maxd));

//    propRanges.push_back(std::make_pair(0,10));
//    propRanges.push_back(std::make_pair(0,65536));
//    propRanges.push_back(std::make_pair(0,65536));

//    propRanges.push_back(std::make_pair(mind,maxd));

    if (p == 0 || p == 3) {
      propRanges.push_back(std::make_pair(mind,maxd));
      propRanges.push_back(std::make_pair(mind,maxd));
    }
}



template<typename I> I static median5(I a, I b, I c, I d, I e)
{
  if (a<b) swap(a,b); // 1
  if (c<d) swap(c,d); // 2
  if (a<c) {swap(a,c); swap(b,d);} // 3
  if (c>e) { // 4
    if (d<e) swap(d,e); // 5
    if (b>d) { // 6
      if (b>c) return c; else return b; // 7
    } else {
      return d; // 6
    } 
  } else {
    if (b>c) { // 5
      if (e>c) { // 6
        if (e>b) return b; else return e; // 7
      } else {
        return c; // 6
      }
    } else {
      return c; // 5
    }
  }
}

#ifdef SMOOTHZOOM
#define imagep(a,b,c,d) (image(a,b,c,d) & PARITYMASK)
#else
#define imagep(a,b,c,d) image(a,b,c,d)
#endif


// Prediction used for interpolation. Does not have to be the same as the guess used for encoding/decoding.
inline ColorVal predict(const Image &image, int z, int p, int r, int c)
{
//    return (p==0?0 : 255);
#ifndef SMOOTHZOOM
#ifdef FIRSTQUANTIZE
    int np = image.numPlanes()/2;
    if (p >= np) {
//      return 4;
      int op = p-np;
//      ColorVal me = image(op,z,r,c);
      ColorVal nA; // most significant bits of neighbor pixels A and B
      ColorVal nB;
      ColorVal lA; // least significant bits of neighbor pixels A and B
      ColorVal lB;
      if (z%2 == 0) { // filling horizontal lines
       nA = image(op,z,r-1,c);
       lA = image(p,z,r-1,c);
       nB = (r+1 < image.rows(z) ? image(op,z,r+1,c) : grey[op]);
       lB = (r+1 < image.rows(z) ? image(p,z,r+1,c) : grey[p]);
      } else { // filling vertical lines
       nA = image(op,z,r,c-1);
       lA = image(p,z,r,c-1);
       nB = (c+1 < image.cols(z) ? image(op,z,r,c+1) : grey[op]);
       lB = (c+1 < image.cols(z) ? image(p,z,r,c+1) : grey[p]);
      }
      ColorVal topleft = (r>0 && c>0 ? image(op,z,r-1,c-1) : grey[op]);
      ColorVal topright = (r>0 && c+1 < image.cols(z) ? image(op,z,r-1,c+1) : grey[op]);
      ColorVal bottomleft = (r+1 < image.rows(z) && c>0 ? image(op,z,r+1,c-1) : grey[op]);
      ColorVal bottomright = (r+1 < image.rows(z) && c+1 < image.cols(z) ? image(op,z,r+1,c+1) : grey[op]);
      topleft <<= QUANTIZATIONSHIFT;
      topright <<= QUANTIZATIONSHIFT;
      bottomleft <<= QUANTIZATIONSHIFT;
      bottomright <<= QUANTIZATIONSHIFT;
      topleft += (r>0 && c>0 ? image(p,z,r-1,c-1) : grey[p]);
      topright += (r>0 && c+1 < image.cols(z) ? image(p,z,r-1,c+1) : grey[p]);
      bottomleft += (r+1 < image.rows(z) && c>0 ? image(p,z,r+1,c-1) : grey[p]);
      bottomright += (r+1 < image.rows(z) && c+1 < image.cols(z) ? image(p,z,r+1,c+1) : grey[p]);
      ColorVal diag1 = topleft+bottomright;
      ColorVal diag2 = topright+bottomleft;

//      if (me >= std::min(nA,nB) && me <= std::max(nA,nB)) {
//      if (nA == nB)
//         return (lA+lB)/2;
//      else
//         return me&QUANTIZATIONMASK;
         nA <<= QUANTIZATIONSHIFT;
         nB <<= QUANTIZATIONSHIFT;
         nA += lA;
         nB += lB;
//         ColorVal avg = (nA + nB)/2;
         ColorVal avg = (nA + nB + diag1 + diag2)/6;
         return (avg&QUANTIZATIONMASK);
//      } else {
//         return (me&QUANTIZATIONMASK);
//      }
    }
#endif

    ColorVal left = (c>0 ? image(p,z,r,c-1) : grey[p]); // always known if filling vertical lines
    ColorVal top = (r>0 ? image(p,z,r-1,c) : grey[p]); // always known if filling horizontal lines
    ColorVal topleft = (r>0 && c>0 ? image(p,z,r-1,c-1) : grey[p]);
    ColorVal topright = (r>0 && c+1 < image.cols(z) ? image(p,z,r-1,c+1) : grey[p]);
    ColorVal bottomleft = (r+1 < image.rows(z) && c>0 ? image(p,z,r+1,c-1) : grey[p]);

//    ColorVal bottomright = (r+1 < image.rows(z) && c+1 < image.cols(z) ? image(p,z,r+1,c+1) : grey[p]);
//    ColorVal diag1 = (topleft+bottomright)/2;
//    ColorVal diag2 = (topright+bottomleft)/2;
    ColorVal gradientTL = left + top - topleft;
//    ColorVal sum = (topleft + topright + bottomleft + bottomright)/4;
    if (z%2 == 0) { // filling horizontal lines
      //    KKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      //    KKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      //    KKKKKKKKKKKKKKKK?
      //    KKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      //
      //    KKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      ColorVal bottom = (r+1 < image.rows(z) ? image(p,z,r+1,c) : top); //grey[p]);
//      ColorVal right = (topright+bottomright)/2;
      ColorVal gradientBL = left + bottom - bottomleft;
      ColorVal avg = (top + bottom)/2;
//      return median3(sum, top, bottom);
//      return median5(gradientTL, gradientBL, avg, diag1, diag2);
      return median3(gradientTL, gradientBL, avg);
//      return median3(gradientTL, left, avg);

    } else { // filling vertical lines
      //   KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      //   KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      //   KKKKKKKKKKKKKKKKK?K K K K K K K
      //   K K K K K K K K K K K K K K K K
      //   K K K K K K K K K K K K K K K K
      //   K K K K K K K K K K K K K K K K
      ColorVal right = (c+1 < image.cols(z) ? image(p,z,r,c+1) : left); //grey[p]);
//      ColorVal bottom = (bottomleft+bottomright)/2;
      ColorVal gradientTR = right + top - topright;
      ColorVal avg = (left + right      )/2;
//      return median3(sum, left, right);
//      return median5(gradientTL, gradientTR, avg, diag1, diag2);
      return median3(gradientTL, gradientTR, avg);
//      return median3(gradientTL, top, avg);
    }
#endif

#ifdef SMOOTHZOOM
  if (z%2 == 0) {
    ColorVal top = imagep(p,z,r-1,c); // actually avg of pixel to guess and pixel just above
    return top;
/*    if (c > 0 && r > 0 && c+1 < image.cols(z) && r+1 < image.rows(z)) {
      ColorVal tb = (top + imagep(p,z,r+1,c)) / 2;
      ColorVal tlbr = (imagep(p,z,r-1,c-1) + imagep(p,z,r+1,c+1)) / 2;
      ColorVal bltr = (imagep(p,z,r+1,c-1) + imagep(p,z,r-1,c+1)) / 2;
      ColorVal bottom = median3(tb,tlbr,bltr); // median of three different estimates for bottom of pixel to guess
      return (top+bottom)/2;
    } else return top;
*/
/*    ColorVal left = (c>0 ? imagep(p,z,r,c-1) : 0);
    ColorVal topleft = (c>0 ? imagep(p,z,r-1,c-1) : 0); // actually avg of left and topleft
    ColorVal gradient = top + left - topleft;
    return gradient;
*/
//    return median3(left,top,gradient);
  } else {
    ColorVal left = imagep(p,z,r,c-1); // actually avg of pixel to guess and pixel just to the left
    return left;

/*    if (c > 0 && r > 0 && c+1 < image.cols(z) && r+1 < image.rows(z)) {
      ColorVal lr = (left + imagep(p,z,r,c+1)) / 2;
      ColorVal tlbr = (imagep(p,z,r-1,c-1) + imagep(p,z,r+1,c+1)) / 2;
      ColorVal bltr = (imagep(p,z,r+1,c-1) + imagep(p,z,r-1,c+1)) / 2;
      ColorVal right = median3(lr,tlbr,bltr); // median of three different estimates for bottom of pixel to guess
      return (left+right)/2;
    } else return left;
*/
/*    ColorVal top = (r>0 ? imagep(p,z,r-1,c) : 0);
    ColorVal topleft = (r>0 ? imagep(p,z,r-1,c-1) : 0);  // actually avg of left and topleft
    ColorVal gradient = top + left - topleft;
    return gradient;
*/
//    return median3(left,top,gradient);
  }
#endif
/*
  if (z%2 == 0) {
    ColorVal top = image(p,z,r-1,c);
    ColorVal bottom = image(p,z,r+1,c);
    if (top < 0) return bottom;
    if (bottom <0) return top;
    return (top+bottom)/2;
  } else {
    ColorVal left = image(p,z,r,c-1);
    ColorVal right = image(p,z,r,c+1);
    if (left < 0) return right;
    if (right < 0) return left;
    return (left+right)/2;
  }
*/
/*
    ColorVal left = image(p,z,r,c-1);
    ColorVal top = image(p,z,r-1,c);
    ColorVal right = image(p,z,r,c+1);
    ColorVal bottom = image(p,z,r+1,c);
    ColorVal topleft = image(p,z,r-1,c-1);
    ColorVal topright = image(p,z,r-1,c+1);
    ColorVal bottomleft = image(p,z,r+1,c-1);
    ColorVal bottomright = image(p,z,r+1,c+1);
    return (  1*(topleft+bottomright)      //  interpolate \  (diagonal, always available)
            + 1*(bottomleft+topright)      //  interpolate /  (diagonal, always available)
            + 2*(z%2==0 ? (top+bottom)     //  interpolate |  (only for horizontal scanlines)
                        : (left+right))    //  interpolate -  (only for vertical scanlines)
           ) / 8;
*/

}
ColorVal predict_and_calcProps(Properties &properties, const ColorRanges *ranges, const Image &image, const int z, const int p, const int r, const int c, ColorVal &min, ColorVal &max) {
    ColorVal guess;
    int which = 0;
    int index = 0;

    if (p != 3) {
    for (int pp = 0; pp < p; pp++) {
//        properties.push_back(imagep(pp,z,r,c));
        properties[index++] = image(pp,z,r,c);
    }
    if (image.numPlanes()>3) properties[index++] = image(3,z,r,c);
    }
//    ranges->minmax(p,properties,min,max);
//    if (min==max) return min;
//    assert(min<max);

/*
    for (int pp = 0; pp < p; pp++) {
        ColorVal oldguess = predict(image,z,pp,r,c);
        if (oldguess > ranges->max(pp)) oldguess = ranges->max(pp);
        if (oldguess < ranges->min(pp)) oldguess = ranges->min(pp);
        properties.push_back(imagep(pp,z,r,c) - oldguess);
    }
*/
    ColorVal left;
    ColorVal top;
    ColorVal topleft = (r>0 && c>0 ? image(p,z,r-1,c-1) : grey[p]);
    ColorVal topright = (r>0 && c+1 < image.cols(z) ? image(p,z,r-1,c+1) : grey[p]);
    ColorVal bottomleft = (r+1 < image.rows(z) && c>0 ? image(p,z,r+1,c-1) : grey[p]);
//    ColorVal bottomright = (r+1 < image.rows(z) && c+1 < image.cols(z) ? image(p,z,r+1,c+1) : grey[p]);

#ifdef FIRSTQUANTIZE
    int np = image.numPlanes()/2;
    if (p >= np) {
      int op = p-np;
      ColorVal me = image(op,z,r,c);
      ColorVal nA; // most significant bits of neighbor pixels A and B
      ColorVal nB;
      ColorVal lA; // least significant bits of neighbor pixels A and B
      ColorVal lB;
      if (z%2 == 0) { // filling horizontal lines
       nA = image(op,z,r-1,c);
       lA = image(p,z,r-1,c);
       nB = (r+1 < image.rows(z) ? image(op,z,r+1,c) : grey[op]);
       lB = (r+1 < image.rows(z) ? image(p,z,r+1,c) : grey[p]);
      } else { // filling vertical lines
       nA = image(op,z,r,c-1);
       lA = image(p,z,r,c-1);
       nB = (c+1 < image.cols(z) ? image(op,z,r,c+1) : grey[op]);
       lB = (c+1 < image.cols(z) ? image(p,z,r,c+1) : grey[p]);
      }
//      if (me >= std::min(nA,nB) && me <= std::max(nA,nB)) {
//      if (nA == nB)
//         return (lA+lB)/2;
//      else
//         return me&QUANTIZATIONMASK;
         if (me >= std::min(nA,nB) && me <= std::max(nA,nB)) which += 1;
         if (nA == nB) which += 2;
         nA <<= QUANTIZATIONSHIFT;
         nB <<= QUANTIZATIONSHIFT;
         nA += lA;
         nB += lB;
         ColorVal avg = (nA + nB)/2;
         guess = (avg&QUANTIZATIONMASK);
    }
#endif
//    ColorVal diag1 = (topleft+bottomright)/2;
//    ColorVal diag2 = (topright+bottomleft)/2;
    if (z%2 == 0) { // filling horizontal lines
      left = (c>0 ? image(p,z,r,c-1) : grey[p]);
      top = image(p,z,r-1,c);
      ColorVal gradientTL = left + top - topleft;
      ColorVal bottom = (r+1 < image.rows(z) ? image(p,z,r+1,c) : top); //grey[p]);
      ColorVal gradientBL = left + bottom - bottomleft;
      ColorVal avg = (top + bottom)/2;
      guess = median3(gradientTL, gradientBL, avg);
//      guess = median5(gradientTL, gradientBL, avg, diag1, diag2);
      ranges->snap(p,properties,min,max,guess);
      if (guess == avg) which = 0;
      else if (guess == gradientTL) which = 1;
      else if (guess == gradientBL) which = 2;
//      else if (guess == diag1) which = 4;
//      else if (guess == diag2) which = 5;
//        guess = median3(avg, diag1, diag2);
//      properties.push_back(top - bottom);
      properties[index++] = top-bottom;

    } else { // filling vertical lines
      left = image(p,z,r,c-1);
      top = (r>0 ? image(p,z,r-1,c) : grey[p]);
      ColorVal gradientTL = left + top - topleft;
      ColorVal right = (c+1 < image.cols(z) ? image(p,z,r,c+1) : left); //grey[p]);
      ColorVal gradientTR = right + top - topright;
      ColorVal avg = (left + right      )/2;
      guess = median3(gradientTL, gradientTR, avg);
//      guess = median5(gradientTL, gradientTR, avg, diag1, diag2);
      ranges->snap(p,properties,min,max,guess);
      if (guess == avg) which = 0;
      else if (guess == gradientTL) which = 1;
      else if (guess == gradientTR) which = 2;
//      else if (guess == diag1) which = 4;
//      else if (guess == diag2) which = 5;
//        guess = median3(avg, diag1, diag2);
//      properties.push_back(left - right);
      properties[index++] = left-right;
    }
    
    //if (guess<min) guess=min;
    //if (guess>max) guess=max;

//    properties.push_back(guess);
//    properties.push_back(which);
    properties[index++]=guess;
    properties[index++]=which;

/*    if (c > 0 && r > 0) { properties.push_back(left - topleft); properties.push_back(topleft - top); }
                 else   { properties.push_back(0); properties.push_back(0); }

    if (c+1 < image.cols(z) && r > 0) properties.push_back(top - topright);
                 else   properties.push_back(0);*/
    if (c > 0 && r > 0) { properties[index++]=left - topleft; properties[index++]=topleft - top; }
                 else   { properties[index++]=0; properties[index++]=0; }

    if (c+1 < image.cols(z) && r > 0) properties[index++]=top - topright;
                 else   properties[index++]=0;

//    if (c > 0 && r > 0 && c+1 < image.cols(z) && r+1 < image.rows(z)) { properties.push_back(topleft-bottomright); properties.push_back(bottomleft-topright); }
//                 else   { properties.push_back(0); properties.push_back(0); }

//    properties.push_back((z%2)*10);    // horizontal or vertical
//    properties.push_back(r * image.zoom_rowpixelsize(z));
//    properties.push_back(c * image.zoom_colpixelsize(z));


   if (p == 0 || p == 3) {
/*   if (r > 1) properties.push_back(image(p,z,r-2,c)-top);    // toptop - top
         else properties.push_back(0);
   if (c > 1) properties.push_back(image(p,z,r,c-2)-left);    // leftleft - left
         else properties.push_back(0);*/
   if (r > 1) properties[index++]=image(p,z,r-2,c)-top;    // toptop - top
         else properties[index++]=0;
   if (c > 1) properties[index++]=image(p,z,r,c-2)-left;    // leftleft - left
         else properties[index++]=0;
   }
    return guess;
}

int plane_zoomlevels(const Image &image, const int beginZL, const int endZL) {
    return image.numPlanes() * (beginZL - endZL + 1);
}

std::pair<int, int> plane_zoomlevel(const Image &image, const int beginZL, const int endZL, int i) {
    assert(i >= 0);
    assert(i < plane_zoomlevels(image, beginZL, endZL));
    // simple order: interleave planes, zoom in
//    int p = i % image.numPlanes();
//    int zl = beginZL - (i / image.numPlanes());

    // more advanced order: give priority to more important plane(s)
    // assumption: plane 0 is Y, plane 1 is I, plane 2 is Q, plane 3 is perhaps alpha, next planes are not important
    const int max_behind[] = {0, 4, 6, 0, 16, 18, 20, 22};
    int np = image.numPlanes();
    if (np>7) {
      // too many planes, do something simple
      int p = i % image.numPlanes();
      int zl = beginZL - (i / image.numPlanes());
      return std::pair<int, int>(p,zl);
    }
    std::vector<int> czl(np);
    for (int &pzl : czl) pzl = beginZL+1;
    int highest_priority_plane = 0;
    if (np >= 4) highest_priority_plane = 3; // alpha first
    int nextp = highest_priority_plane;
    while (i >= 0) {
      czl[nextp]--;
      i--;
      if (i<0) break;
      nextp=highest_priority_plane;
      for (int p=0; p<np; p++) {
        if (czl[p] > czl[highest_priority_plane] + max_behind[p]) {
//          printf("czl[%i] > %i + %i\n",czl[p] , czl[0] , max_behind[p]);
          nextp = p; break;
        }
      }
      // ensure that nextp is not at the most detailed zoomlevel yet
      while (czl[nextp] <= endZL) nextp = (nextp+1)%np;
    }
    int p = nextp;
    int zl = czl[p];

    return std::pair<int, int>(p,zl);
}

template<typename Coder, typename ParityCoder> void encode_jif2_inner(std::vector<Coder*> &coders, ParityCoder &parityCoder, const Image &image, const ColorRanges *ranges, const int beginZL, const int endZL)
{
    long fs = ftell(f);
    ColorVal min,max;
    int nump = image.numPlanes();
//    for (int p = 0; p < image.numPlanes(); p++) {
//    for (int z = beginZL; z >= endZL; z--) {
#ifdef FIRSTQUANTIZE
    int np = image.numPlanes()/2;
#endif
    for (int i = 0; i < plane_zoomlevels(image, beginZL, endZL); i++) {
      std::pair<int, int> pzl = plane_zoomlevel(image, beginZL, endZL, i);
      int p = pzl.first;
      int z = pzl.second;
      if (ranges->min(p) < ranges->max(p))
        fprintf(stdout,"[%i/%i] ENC(p:%i,z:%i), %ix%i\t",i,plane_zoomlevels(image, beginZL, endZL)-1,p,z,image.rows(z),image.cols(z));
      else continue;
      fflush(stdout);
      Properties properties((nump>3?NB_PROPERTIESA[p]:NB_PROPERTIES[p]));
      if (z % 2 == 0) {
        // horizontal: scan the odd rows, output pixel values
          for (int r = 1; r < image.rows(z); r += 2) {
            for (int c = 0; c < image.cols(z); c++) {
//                    properties.clear();
                    if (nump>3 && p<3 && image(3,z,r,c) == 0) continue;
                    ColorVal guess = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
                    ColorVal curr = image(p,z,r,c);
                    //ranges->position(p,properties,image(p,z,r,c));
                    assert (curr <= max);
                    assert (curr >= min);
                    assert (ranges->reconstruct(p,properties,curr) == image(p,z,r,c));
#ifdef SMOOTHZOOM
                    int parity = (curr & PARITYBIT ? 1 : 0);
                    curr &= PARITYMASK;
                    parityCoder.write(parity);
#endif
#ifdef FIRSTQUANTIZE
                    if (p>=np) coders[p]->write_int(properties, QUANTIZATIONSHIFT, curr);
                    else
#endif
                    coders[p]->write_int(properties, min - guess, max - guess, curr - guess);
            }
          }
      } else {
        // vertical: scan the odd columns
          for (int r = 0; r < image.rows(z); r++) {
            for (int c = 1; c < image.cols(z); c += 2) {
//                    properties.clear();
                    if (nump>3 && p<3 && image(3,z,r,c) == 0) continue;
                    ColorVal guess = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
                    ColorVal curr = image(p,z,r,c);
                    assert (curr <= max);
                    assert (curr >= min);
                    assert (ranges->reconstruct(p,properties,curr) == image(p,z,r,c));
#ifdef SMOOTHZOOM
                    int parity = (curr & PARITYBIT ? 1 : 0);
                    curr &= PARITYMASK;
                    parityCoder.write(parity);
#endif
#ifdef FIRSTQUANTIZE
                    if (p>=np) coders[p]->write_int(properties, QUANTIZATIONSHIFT, curr);
                    else
#endif
                    coders[p]->write_int(properties, min - guess, max - guess, curr - guess);
            }
          }
      }
      long nfs = ftell(f);
      long pixels = image.cols(z)*image.rows(z)/2;
      if (nfs-fs > 0) fprintf(stdout,"filesize:%li (+%li bytes, %li pixels, %f bpp)\n", nfs, nfs-fs, pixels, 8.0*(nfs-fs)/pixels );
      fs = nfs;
    }
}

template<typename Rac, typename Coder> void encode_jif2_pass(Rac &rac, const Image &image, const ColorRanges *ranges, std::vector<Tree> &forest, const int beginZL, const int endZL, int repeats)
{
    std::vector<Coder*> coders;
    for (int p = 0; p < ranges->numPlanes(); p++) {
        Ranges propRanges;
        initPropRanges(propRanges, *ranges, p);
        coders.push_back(new Coder(rac, propRanges, forest[p]));
    }

    if (beginZL == image.zooms()) {
    // special case: very left top pixel must be written first to get it all started
    SimpleSymbolCoder<JifBitChanceMeta, Rac, 24> metaCoder(rac);
    // in case of smooth zoom, this pixel is the average of the entire image
#ifdef SMOOTHZOOM
    printf("Summary pixel: ");
#else
    printf("First pixel: ");
#endif
    for (int p = 0; p < image.numPlanes(); p++) {
        ColorVal curr = image(p,0,0);
        metaCoder.write_int(ranges->min(p), ranges->max(p), curr);
        printf("p%i=%i, ",p,image(p,0,0));
    }
    printf("\n");
    }
    SimpleBitCoder<JifBitChanceParities,Rac> parityCoder(rac);

    if (repeats>1) printf("Iterating %i times to find a better tree.\n",repeats);
    while(repeats-- > 0) {
     encode_jif2_inner(coders, parityCoder, image, ranges, beginZL, endZL);
    }

    for (int p = 0; p < image.numPlanes(); p++) {
        coders[p]->simplify();
    }

    for (int p = 0; p < image.numPlanes(); p++) {
#ifdef STATS
        indent(0); printf("Plane %i\n", p);
        coders[p]->info(0+1);
#endif
        delete coders[p];
    }
}

void encode_jif2_interpol_zero_alpha(Image &image, const ColorRanges *ranges, const int beginZL, const int endZL)
{
    ColorVal min,max;
    printf("Replacing fully transparent pixels with predicted pixel values at the other planes\n");
    for (int i = 0; i < plane_zoomlevels(image, beginZL, endZL); i++) {
      std::pair<int, int> pzl = plane_zoomlevel(image, beginZL, endZL, i);
      int p = pzl.first;
      int z = pzl.second;
      if (p == 3) continue;
      Properties properties(NB_PROPERTIES[p]);
      if (z % 2 == 0) {
        // horizontal: scan the odd rows
          for (int r = 1; r < image.rows(z); r += 2) {
            for (int c = 0; c < image.cols(z); c++) {
               if (image(3,z,r,c) == 0) image(p,z,r,c) = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
            }
          }
      } else {
        // vertical: scan the odd columns
          for (int r = 0; r < image.rows(z); r++) {
            for (int c = 1; c < image.cols(z); c += 2) {
               if (image(3,z,r,c) == 0) image(p,z,r,c) = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
            }
          }
      }
    }
}


// interpolate rest of the image
// used when decoding lossy
void decode_jif2_inner_interpol(Image &image, const ColorRanges *ranges, int I, const int beginZL, const int endZL, int R)
{
//    Properties properties;
//    ColorVal min,max;
    for (int i = I; i < plane_zoomlevels(image, beginZL, endZL); i++) {
      std::pair<int, int> pzl = plane_zoomlevel(image, beginZL, endZL, i);
      int p = pzl.first;
      int z = pzl.second;
      fprintf(stdout,"[%i/%i] INTERPOLATE Plane %i, Zoomlevel %i, size: %ix%i\n",i,plane_zoomlevels(image, beginZL, endZL)-1,p,z,image.rows(z),image.cols(z));
      if (z % 2 == 0) {
        // horizontal: scan the odd rows
          for (int r = (I==i?R:1); r < image.rows(z); r += 2) {
            for (int c = 0; c < image.cols(z); c++) {
               image(p,z,r,c) = predict(image,z,p,r,c);
//                 image(p,z,r,c) = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
//               image(p,z,r,c) = ranges->snap(p,predict(image,z,p,r,c),z,r,c);
                // snapping can cause artifacts, better snap at the end
            }
          }
      } else {
        // vertical: scan the odd columns
          for (int r = (I==i?R:0); r < image.rows(z); r++) {
            for (int c = 1; c < image.cols(z); c += 2) {
               image(p,z,r,c) = predict(image,z,p,r,c);
//                 image(p,z,r,c) = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
//               image(p,z,r,c) = ranges->snap(p,predict(image,z,p,r,c),z,r,c);
            }
          }
      }
    }
}

template<typename Coder, typename ParityCoder> void decode_jif2_inner(std::vector<Coder*> &coders, ParityCoder &parityCoder, Image &image, const ColorRanges *ranges, const int beginZL, const int endZL, int lastI)
{
    ColorVal min,max;
    int nump = image.numPlanes();
#ifdef FIRSTQUANTIZE
    int np = image.numPlanes()/2;
#endif
    // decode
    for (int i = 0; i < plane_zoomlevels(image, beginZL, endZL); i++) {
      std::pair<int, int> pzl = plane_zoomlevel(image, beginZL, endZL, i);
      int p = pzl.first;
      int z = pzl.second;
      if (lastI != -1  && i > lastI) {
              decode_jif2_inner_interpol(image, ranges, i, beginZL, endZL, (z%2 == 0 ?1:0));
              return;
      }
      if (ranges->min(p) < ranges->max(p))
        fprintf(stdout,"[%i/%i] DEC Plane %i, Zoomlevel %i, size: %ix%i\n",i,plane_zoomlevels(image, beginZL, endZL)-1,p,z,image.rows(z),image.cols(z));
      else continue;
      ColorVal curr;
      Properties properties((nump>3?NB_PROPERTIESA[p]:NB_PROPERTIES[p]));
      if (z % 2 == 0) {
          for (int r = 1; r < image.rows(z); r += 2) {
#ifdef CHECK_FOR_BROKENFILES
            if (feof(f)) {
              printf("Unexpected file end. Interpolation from now on.\n");
              decode_jif2_inner_interpol(image, ranges, i, beginZL, endZL, r);
              return;
            }
#endif
            for (int c = 0; c < image.cols(z); c++) {
                     ColorVal guess = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
                     if (nump>3 && p<3 && image(3,z,r,c) == 0) {image(p,z,r,c) = guess; continue;}
#ifdef SMOOTHZOOM
                     int parity = parityCoder.read();
#endif
#ifdef FIRSTQUANTIZE
                    if (p>=np) curr = coders[p]->read_int(properties, QUANTIZATIONSHIFT);
                    else
#endif
                     curr = coders[p]->read_int(properties, min - guess, max - guess) + guess;
                     image(p,z,r,c) = curr; //ranges->reconstruct(p,properties,curr);
#ifdef SMOOTHZOOM
                     if (parity) image(p,z,r,c) += PARITYBIT;
#endif
            }
        }
      } else {
          for (int r = 0; r < image.rows(z); r++) {
#ifdef CHECK_FOR_BROKENFILES
            if (feof(f)) {
              printf("Unexpected file end. Interpolation from now on.\n");
              decode_jif2_inner_interpol(image, ranges, i, beginZL, endZL, r);
              return;
            }
#endif
            for (int c = 1; c < image.cols(z); c += 2) {
                     ColorVal guess = predict_and_calcProps(properties,ranges,image,z,p,r,c,min,max);
                     if (nump>3 && p<3 && image(3,z,r,c) == 0) {image(p,z,r,c) = guess; continue;}
#ifdef SMOOTHZOOM
                     int parity = parityCoder.read();
#endif
#ifdef FIRSTQUANTIZE
                    if (p>=np) curr = coders[p]->read_int(properties, QUANTIZATIONSHIFT);
                    else
#endif
                     curr = coders[p]->read_int(properties, min - guess, max - guess) + guess;
                     image(p,z,r,c) = curr; //ranges->reconstruct(p,properties,curr);
#ifdef SMOOTHZOOM
                     if (parity) image(p,z,r,c) += PARITYBIT;
#endif
            }
        }
      }
    }
}

template<typename Rac, typename Coder> void decode_jif2_pass(Rac &rac, Image &image, const ColorRanges *ranges, std::vector<Tree> &forest, const int beginZL, const int endZL, int lastI)
{
    std::vector<Coder*> coders;
    for (int p = 0; p < image.numPlanes(); p++) {
        Ranges propRanges;
        initPropRanges(propRanges, *ranges, p);
        coders.push_back(new Coder(rac, propRanges, forest[p]));
    }

    if (beginZL == image.zooms()) {
    // special case: very left top pixel must be read first to get it all started
    SimpleSymbolCoder<JifBitChanceMeta, Rac, 24> metaCoder(rac);
#ifdef SMOOTHZOOM
    printf("Summary pixel: ");
#else
    printf("First pixel: ");
#endif
    for (int p = 0; p < image.numPlanes(); p++) {
        image(p,0,0) = metaCoder.read_int(ranges->min(p), ranges->max(p));
        printf("p%i=%i, ",p,image(p,0,0));
    }
    printf("\n");
    }

    SimpleBitCoder<JifBitChanceParities,Rac> parityCoder(rac);

    decode_jif2_inner(coders, parityCoder, image, ranges, beginZL, endZL, lastI);

    for (int p = 0; p < image.numPlanes(); p++) {
        delete coders[p];
    }
}


/******************************************/
/*   General encoding/decoding            */
/******************************************/

template<typename BitChance, typename Rac> void encode_tree(Rac &rac, const ColorRanges *ranges, const std::vector<Tree> &forest, const int encoding)
{
    for (int p = 0; p < ranges->numPlanes(); p++) {
        Ranges propRanges;
        if (encoding==1) initPropRanges_ffv1(propRanges, *ranges, p);
        else initPropRanges(propRanges, *ranges, p);
        MetaPropertySymbolCoder<BitChance, Rac> metacoder(rac, propRanges);
//        forest[p].print(stdout);
        metacoder.write_tree(forest[p]);
    }
}
template<typename BitChance, typename Rac> void decode_tree(Rac &rac, const ColorRanges *ranges, std::vector<Tree> &forest, const int encoding)
{
    for (int p = 0; p < ranges->numPlanes(); p++) {
        Ranges propRanges;
        if (encoding==1) initPropRanges_ffv1(propRanges, *ranges, p);
        else initPropRanges(propRanges, *ranges, p);
        MetaPropertySymbolCoder<BitChance, Rac> metacoder(rac, propRanges);
        metacoder.read_tree(forest[p]);
//        forest[p].print(stdout);
    }
}


bool encode(const char* filename, Image &image, std::vector<std::string> transDesc, int encoding)
{
// not computing checksum until after transformations and potential zero-alpha changes
//    uint32_t checksum = image.checksum();


    f = fopen(filename,"w");
    RacOut rac(f);

    switch(encoding) {
        case 1: write_name(rac, "JIF1"); break;
        case 2: write_name(rac, "JIF2"); break;
        default: fprintf(stderr,"Unknown encoding: %i\n", encoding); return false;
    }

    SimpleSymbolCoder<JifBitChanceMeta, RacOut, 24> metaCoder(rac);
    int numPlanes = image.numPlanes();
    metaCoder.write_int(1, 16, numPlanes);
    metaCoder.write_int(1, 65536, image.cols());
    metaCoder.write_int(1, 65536, image.rows());
    printf("Input planes: ");
    for (int p = 0; p < numPlanes; p++) {
//        metaCoder.write_int(1, 64, image.subSampleR(p));
//        metaCoder.write_int(1, 64, image.subSampleC(p));
//        metaCoder.write_int(-16777216, 16777215, image.min(p));
//        metaCoder.write_int(0, 16777216, image.max(p) - image.min(p));
        assert(image.min(p) == 0);
        metaCoder.write_int(1, 16, ilog2(image.max(p)+1));
        printf("[%i] %i bpp (%i..%i) \t",p,ilog2(image.max(p)+1),image.min(p), image.max(p));
    }
    printf("\n");

    std::vector<const ColorRanges*> rangesList;
    std::vector<Transform*> transforms;
    rangesList.push_back(getRanges(image));
    for (unsigned int i=0; i<transDesc.size(); i++) {
        Transform *trans = create_transform(transDesc[i]);
        if (!trans->init(rangesList.back()) || !trans->process(rangesList.back(), image)) {
            fprintf(stderr, "Transform '%s' failed\n", transDesc[i].c_str());
        } else {
            printf("Doing transform '%s'\n", transDesc[i].c_str());
            rac.write(true);
            write_name(rac, transDesc[i]);
            trans->save(rangesList.back(), rac);
            rangesList.push_back(trans->meta(image, rangesList.back()));
            trans->data(image);
        }
    }
    rac.write(false);
//    const ColorRanges* ranges = rangesList[1];
    const ColorRanges* ranges = rangesList.back();
    grey.clear();
    for (int p = 0; p < ranges->numPlanes(); p++) grey.push_back((ranges->min(p)+ranges->max(p))/2);

    int mbits = 0;
    for (int p = 0; p < ranges->numPlanes(); p++) {
        int nBits = ilog2((ranges->max(p) - ranges->min(p))*2-1)+1;
        if (nBits > mbits) mbits = nBits;
    }
    const int bits = 10;
    if (mbits > bits) { printf("OOPS: %i > %i\n",mbits,bits); return false;}

    // two passes
    std::vector<Tree> forest(ranges->numPlanes(), Tree());
    RacDummy dummy;

    int roughZL = 0;
    if (ranges->numPlanes() > 3) switch(encoding) {
        case 1: encode_ffv1_interpol_zero_alpha(image, ranges); break;
        case 2: encode_jif2_interpol_zero_alpha(image, ranges, image.zooms(), 0); break;
    }

    uint32_t checksum = image.checksum();

    if (encoding == 2) {
      roughZL = image.zooms() - NB_NOLEARN_ZOOMS;
      if (roughZL < 0) roughZL = 0;
      fprintf(stdout,"Encoding rough data\n");
      encode_jif2_pass<RacOut, FinalPropertySymbolCoder<JifBitChancePass2, RacOut, bits> >(rac, image, ranges, forest, image.zooms(), roughZL+1, 1);
    }

    fprintf(stdout,"Encoding data (pass 1)\n");
    switch(encoding) {
        case 1: encode_ffv1_pass<RacDummy, PropertySymbolCoder<JifBitChancePass1, RacDummy, bits> >(dummy, image, ranges, forest, TREE_LEARN_REPEATS); break;
        case 2: encode_jif2_pass<RacDummy, PropertySymbolCoder<JifBitChancePass1, RacDummy, bits> >(dummy, image, ranges, forest, roughZL, 0, TREE_LEARN_REPEATS); break;
    }

    fprintf(stdout,"Encoding tree\n");
    long fs = ftell(f);
    encode_tree<JifBitChanceTree, RacOut>(rac, ranges, forest, encoding);
    fprintf(stdout,"Rough data total: %li bytes.  Tree total: %li bytes.\n", fs, ftell(f)-fs);
    fprintf(stdout,"Encoding data (pass 2)\n");
    switch(encoding) {
        case 1: encode_ffv1_pass<RacOut, FinalPropertySymbolCoder<JifBitChancePass2, RacOut, bits> >(rac, image, ranges, forest, 1); break;
        case 2: encode_jif2_pass<RacOut, FinalPropertySymbolCoder<JifBitChancePass2, RacOut, bits> >(rac, image, ranges, forest, roughZL, 0, 1); break;
    }

/*
    if (encoding == 2) {
     std::vector<Tree> forest2(ranges->numPlanes(), Tree());
     fprintf(stdout,"Encoding data (pass 3)\n");
     encode_jif2_pass<RacDummy, PropertySymbolCoder<JifBitChancePass1, RacDummy> >(dummy, image, ranges, forest2, 1, 0);
     fprintf(stdout,"Encoding tree 2\n");
     long fs = ftell(f);
     encode_tree<JifBitChanceTree, RacOut>(rac, ranges, forest2);
     fprintf(stdout,"Medium data total: %li bytes.  Tree 2 total: %li bytes.\n", fs, ftell(f)-fs);
     fprintf(stdout,"Encoding data (pass 4)\n");
     encode_jif2_pass<RacOut, FinalPropertySymbolCoder<JifBitChancePass2, RacOut> >(rac, image, ranges, forest2, 1, 0);
    }
*/
    fprintf(stdout,"Encoding done\n");
//    rac.flush();
    fprintf(stdout,"Writing checksum: %X\n", checksum);
    metaCoder.write_int(0, 0xFFFF, checksum / 0x10000);
    metaCoder.write_int(0, 0xFFFF, checksum & 0xFFFF);
    rac.flush();
    fclose(f);

//    fprintf(stdout,"Cleaning up\n");
    for (int i=transforms.size()-1; i>=0; i--) {
        delete transforms[i];
    }
    transforms.clear();
    for (unsigned int i=0; i<rangesList.size(); i++) {
        delete rangesList[i];
    }
    rangesList.clear();

    return true;
}



bool decode(const char* filename, Image &image, int lastI)
{
    image.reset();

    f = fopen(filename,"r");
    RacIn rac(f);
    int encoding=0;

    std::string str = read_name(rac);
    if (str == "JIF1") {
        encoding=1;
    } else if (str == "JIF2") {
        encoding=2;
    } else {
        fprintf(stderr,"Unknown magic '%s'\n", str.c_str());
        return false;
    }

    SimpleSymbolCoder<JifBitChanceMeta, RacIn, 24> metaCoder(rac);
    int numPlanes = metaCoder.read_int(1, 16);
    int width = metaCoder.read_int(1, 65536);
    int height = metaCoder.read_int(1, 65536);
    image.init(width, height, 0, 0, 0);
    for (int p = 0; p < numPlanes; p++) {
        int subSampleR = 1; //metaCoder.read_int(1, 64);
        int subSampleC = 1; //metaCoder.read_int(1, 64);
//        int min = metaCoder.read_int(-16777216, 16777215);
//        int max = metaCoder.read_int(0, 16777216) + min;
        int min = 0;
        int max = (1 << metaCoder.read_int(1, 16)) - 1;
        image.add_plane(min, max, subSampleR, subSampleC);
        printf("plane %i: %i bits per pixel (%i..%i)\n",p,ilog2(image.max(p)+1),image.min(p), image.max(p));
    }

    std::vector<const ColorRanges*> rangesList;
    std::vector<Transform*> transforms;
    rangesList.push_back(getRanges(image));
    while (rac.read()) {
        std::string desc = read_name(rac);
        Transform *trans = create_transform(desc);
        if (!trans) {
            fprintf(stderr,"Unknown transformation '%s'\n", desc.c_str());
            return false;
        }
        if (!trans->init(rangesList.back())) {
            fprintf(stderr,"Transformation '%s' failed\n", desc.c_str());
            return false;
        }
// todo: fix this (if smooth zooming would be used again, it's not a good idea for lossless compression though!)
//        if (desc == "ZOOM") trans->configure(zoomlevel);
        printf("Doing transform '%s'\n", desc.c_str());
        trans->load(rangesList.back(), rac);
        rangesList.push_back(trans->meta(image, rangesList.back()));
        transforms.push_back(trans);
    }
//    const ColorRanges* ranges = rangesList[1];
    const ColorRanges* ranges = rangesList.back();
    grey.clear();
    for (int p = 0; p < ranges->numPlanes(); p++) grey.push_back((ranges->min(p)+ranges->max(p))/2);

    for (int p = 0; p < numPlanes; p++) {
        if (ranges->min(p) >= ranges->max(p)) {
             printf("Constant plane %i at color value %i\n",p,ranges->min(p));
             for (ColorVal& x : image(p).data) x=ranges->min(p);
        }
    }

    int mbits = 0;
    for (int p = 0; p < ranges->numPlanes(); p++) {
        int nBits = ilog2((ranges->max(p) - ranges->min(p))*2-1)+1;
        if (nBits > mbits) mbits = nBits;
    }
    const int bits = 10;
    if (mbits > bits) { printf("OOPS: %i > %i\n",mbits,bits); return false;}


    std::vector<Tree> forest(ranges->numPlanes(), Tree());

    int roughZL = 0;
    if (encoding == 2) {
      roughZL = image.zooms() - NB_NOLEARN_ZOOMS;
      if (roughZL < 0) roughZL = 0;
      fprintf(stdout,"Decoding rough data\n");
      decode_jif2_pass<RacIn, FinalPropertySymbolCoder<JifBitChancePass2, RacIn, bits> >(rac, image, ranges, forest, image.zooms(), roughZL+1, -1);
    }
    if (encoding == 2 && lastI < -2) {
      printf("Not decoding tree\n");
    } else {
      fprintf(stdout,"Decoding tree\n");
      decode_tree<JifBitChanceTree, RacIn>(rac, ranges, forest, encoding);
    }
    switch(encoding) {
        case 1: fprintf(stdout,"Decoding data (FFV1)\n");
                decode_ffv1_pass<RacIn, FinalPropertySymbolCoder<JifBitChancePass2, RacIn, bits> >(rac, image, ranges, forest);
                break;
        case 2: fprintf(stdout,"Decoding data (JIF2)\n");
                decode_jif2_pass<RacIn, FinalPropertySymbolCoder<JifBitChancePass2, RacIn, bits> >(rac, image, ranges, forest, roughZL, 0, lastI);
                break;
    }
    fprintf(stdout,"Decoding done, read %ld bytes\n",ftell(f));


    if (lastI < 0) {
      uint32_t checksum = image.checksum();
      fprintf(stdout,"Computed checksum: %X\n", checksum);
      uint32_t checksum2 = metaCoder.read_int(0, 0xFFFF);
      checksum2 *= 0x10000;
      checksum2 += metaCoder.read_int(0, 0xFFFF);
      fprintf(stdout,"Read checksum: %X\n", checksum2);
      if (checksum != checksum2) printf("\nCORRUPTION DETECTED!\n");
    } else {
      fprintf(stdout,"Not checking checksum, lossy partial decoding was chosen.\n");
    }

    for (int i=transforms.size()-1; i>=0; i--) {
        transforms[i]->invData(image);
        delete transforms[i];
    }
    transforms.clear();


    for (unsigned int i=0; i<rangesList.size(); i++) {
        delete rangesList[i];
    }
    rangesList.clear();

    fclose(f);
    return true;
}


int main(int argc, char **argv)
{
    Image image;
    if (argc == 3) {
        image.load(argv[1]);
        std::vector<std::string> desc;
//        desc.push_back("SGR");
        desc.push_back("YIQ");
//        desc.push_back("BSH");
#ifdef FIRSTQUANTIZE
        desc.push_back("QTZ");
#endif
#ifdef SMOOTHZOOM
        desc.push_back("ZOOM");
#endif
#ifdef PERMUTEPLANES
        desc.push_back("PLP");
#endif
        desc.push_back("BND");
        desc.push_back("ACB");
        encode(argv[2], image, desc, 1);
    } else if (argc == 4) {
        decode(argv[2], image, -1);
        image.save(argv[3]);
    } else if (argc == 5) {
        decode(argv[3], image, strtol(argv[2],NULL,0));
        image.save(argv[4]);
    } else {
        fprintf(stderr, "Usage: %s [-d [level]] source dest\n", argv[0]);
        return 1;
    }
    return 0;
}
