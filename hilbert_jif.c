/*  jif - an experimental segment-based image format encoder/decoder
 Copyright (C) 2010  Jon Sneyers & Pieter Wuille

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 Parts of this code are based on code from the FFMPEG project, in
 particular:
 - ffv1.c - Copyright (c) 2003 Michael Niedermayer <michaelni@gmx.at>
 - common.h - copyright (c) 2006 Michael Niedermayer <michaelni@gmx.at>
 - rangecoder.c - Copyright (c) 2004 Michael Niedermayer <michaelni@gmx.at>
 */

/*
 TODO:
 =====

 * Asymmetrische splits / smalle splits
 * Parent diff als context voor child diffs
 * Andere interpolaties
   - bicubisch?
   - specifiek voor lossless: op grote schaal geen dure interpolaties proberen
 * Encoden in twee fases: 
        - eerst beste interpolaties zoeken (encoden met alle mogelijke interpolaties en counts bijhouden)
        - dan "echt" encoderen waarbij interpolatie-methoden die minder dan X % gebruikt worden weggelaten worden
        - lijst van methodes in file header (volgorde belangrijk!)
 * Andere kleurconfiguraties (YUV, greyscale, RGBA, 16-bit)
 * Subsampling
 * Breadth-first
 * beslissen welke range encoder (16-bit, 24-bit, 40-bit of 56-bit)
 * accurate fraction encoder optimaliseren
 * Estimator verbeteren (mediaan?)
 * Area-distance functie
 * Meta-sectie?
   * EXIF tags includen?
   * Comments
   * Zlib gecomprimeerd, Rac gecomprimeerd, ongecomprimeerd?
 * Betere header?
   * Afbeeldingsgrootte
   * Kleurschema (YIQ/YUV/RGB)
   * Subsampling
   * Kleurdiepte van kanalen
   * Interpolatie-algoritme
   * Predictie-algoritme
   * Aanwezigheid van CRC, of CRC zelf al
   * Grootte van RAC chunk (nu is het gewoon 'tot einde van file')
   * Uitbreidbaarheid in gedachte houden...
 * Naam verzinnen voor dit formaat
   -> "RRS" ?  (Recursive Rectangle Splitting)
   -> "RRR" ?  (Recursive Rectangles with Range encoding)
   -> "SIR" ?  (Split, Interpolate, Recurse)
   -> "SIRF" ? (Splitting, Interpolating, Recursing Format)
   -> "SIRI" ? (Split, Improve, Recurse, Interpolate)
   -> "RIF" ?  (Recursively Interpolating Format)
   -> "RIS" ?  (Recurisve Image Subdivision)
   -> "JPEG ?  (Jon & Pieter's Encoding of Graphics)
   -> "SWF" ?  (Sneyers-Wuille Format)

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "config.h"
#include "util.h"
#include "symbol.h"
#include "image/pixel.h"
#include "image/image.h"
#include "color.h"
#include "interpol.h"
#include "crc32k.h"
#include "log4k.h"

// magic header bytes
uint8_t const magic[] = "JiF1";

// magic marker at the end (inside rac)
uint16_t const imagic = 0x5E2D;

// full context table for range encoder
typedef struct {
  chs_t splitCont[QZ_SIZE_S][3][4][QZ_BVAR_S]; // 0=split and continue      1=stop & interpolate
//  chs_t splitContA[QZ_SIZE_S][4][QZ_BVAR_S]; // 0=split and continue      1=stop & interpolate
  chs_t splitShort[32]; // 1=split along short axis  0=split along long axis and stop
  chs_t splitAssym[32]; // 1=split assymetrically    0=split centrally
  chs_t splitSided[32]; // 1=split near high coords  0=split near low coords

  chs_t interpolationMethod[32][3][MAX_INTERPOLATIONS][MAX_INTERPOLATIONS]; 

  // symbol chances for each size subpicture and color plane
  symb_chs_t diff[QZ_SIZE][3][QZ_LRDIFF][QZ_GRAD][QZ_CROSS][QZ_VAR][QZ_BVAR];
//  symb_chs_t mask[8][QZ_SIZE_S][QZ_BVAR_S][QZ_BVAR_S][QZ_BVAR_S];

  symb_chs_t methods;
  symb_chs_t r_min[3];
  symb_chs_t r_max[3];

  uint64_t bitsSplit[4]; uint32_t symbSplit[4];
  uint64_t bitsInterpol[3]; uint32_t symbInterpol[3];
  uint64_t bitsPixeldata[3]; uint32_t symbPixeldata[3];
} ctx_t;

// data for encoder/decoder
typedef struct {
  image_t *source; // for encoder: source image, for decoder: decoded image
  image_t *rec; // image for reconstruction
  int maxD[3]; // maximum color distances (encoder only)
  int qz[3]; // quantization (1 for lossless)
  double factor; // allowed factor (encoder only)
  symb_coder_t coder; // symbol encoder private data
  int range_min[3];  //
  int range_max[3];  // ranges for Y,I,Q

  ctx_t *ctx; // context table
  
  int nb_interpolation_methods;
  int interpolation_methods[MAX_INTERPOLATIONS];

  // statistics
  int outputted_pixels[3]; // number of outputted pixels
  int splits[3]; // number of splits in each color plane
  int pLeft,pRight,pTop,pBottom;
  int method[MAX_INTERPOLATIONS+1];
  uint64_t method_gain[MAX_INTERPOLATIONS+1];
  int minsplit;
} coder_t;

// initial color
static const pixel_t purple = {.d= {240,500,0},.a=0};

// statistics
double static maxdist[3] = {0,0,0};
long int static stats_bvar[3][QZ_BVAR] = {};

// plane names
static const char planesRGB[4] = "RGBA";
static const char planesYIQ[4] = "YIQA";

uint16_t path[3] = {STRETCH(0),STRETCH(255),STRETCH(255)};
int pixelcounter = 0;

int static inline clamp(int val, int min, int max) {
  if (val<min) val=min;
  if (val>max) val=max;
  return val;
}

// division rounding down (instead of toward zero)
int static inline rdiv(int num, int den) {
  if (num>=0) {
    return num/den;
  } else {
    return (num-den+1)/den;
  }
}

//#include "init_chances.h"

int static inline quantize_log(int val, int maxval, int sign, int levels) {
  if (levels==1 || val==0) return 0;
  int rlevels=(sign ? (levels-1)/2 : levels-1);
  int pre=(log4k[abs(4096/val)]*rlevels)/(log4k[4096/maxval-1]);
  int result = 1+pre;
  if (sign) {
    return 1+pre*2+(val<0);
  } else {
//    fprintf(stderr,"val:%i, maxval:%i, levels:%i, pre:%i\n",val,maxval,levels,pre);
    if (result >= levels) result=levels-1;
    assert(result < levels);
    assert(result > 0);
    return result;
  }
}

int static inline quantize_log_uint32(uint32_t val, int levels) {
  if (val==0 || levels==1) return 0;
  levels--;
  int add=0;
  while (val < (1UL<<31)) {
    val *= 2;
    add++;
  }
  int app=(val >= 0xFFF80000) ? 4096 : (val+(1<<19))/(1<<20);
  int l=log4k[app]+add*log4k[2048]; // 0..174752 (5461 voor val=1<<31)
  int pre=((l-1)*levels)/174752;
  return levels-pre;
}
int qbvar(uint32_t var, int quant) {
    int result=quantize_log_uint32(2*var+10,quant+4)-4;
    if (result<0) result=0;
    return result;
}

void static ctx_init(ctx_t *ctx, int cutoff) {

  // initialize split chances
  ctx->bitsSplit[3]=0;
  ctx->symbSplit[3]=0;
  for (int c = 0; c < 3; c++) {
    ctx->bitsSplit[c]=0;
    ctx->symbSplit[c]=0;
    ctx->bitsInterpol[c]=0;
    ctx->symbInterpol[c]=0;
    ctx->bitsPixeldata[c]=0;
    ctx->symbPixeldata[c]=0;
    for (int i = 0; i < QZ_SIZE_S; i++) {
      for (int j = 0; j < 4; j++) {
        for (int n = 0; n < QZ_BVAR_S; n++) {
//          chs_init(&ctx->splitCont[i][c][j][n],clamp(96+i*(25)*(c+3)+100*n/QZ_BVAR_S,cutoff,4096-cutoff));
          chs_init(&ctx->splitCont[i][c][j][n],clamp(50*(c+5)*(QZ_SIZE_S-i)/QZ_SIZE_S + 25*(c) + 50*(c+5)*(QZ_BVAR_S-n)/QZ_BVAR_S,cutoff,4096-cutoff));
//          chs_init(&ctx->splitCont[i][c][j][n],2048);
        }
      }
      for (int j = 0; j < MAX_INTERPOLATIONS; j++) {
        for (int k = 0; k < MAX_INTERPOLATIONS; k++) {
          chs_init(&ctx->interpolationMethod[i][c][j][k],2048);
        }
      }
    }
  }

/*
    for (int i = 0; i < QZ_SIZE_S; i++) {
      for (int j = 0; j < 4; j++) {
        for (int n = 0; n < QZ_BVAR_S; n++) {
//          chs_init(&ctx->splitContA[i][j][n],clamp(1000*(QZ_SIZE_S-i)/QZ_SIZE_S + 1000*(QZ_BVAR_S-n)/QZ_BVAR_S,cutoff,4096-cutoff));
          chs_init(&ctx->splitContA[i][j][n],2048);
        }
      }
    }
*/

  // initialize symbol chances
  for (int i = 0; i < QZ_SIZE; i++) 
  for (int c = 0; c < 3; c++) 
  for (int j = 0; j < QZ_LRDIFF; j++) 
  for (int k = 0; k < QZ_GRAD; k++) 
  for (int l = 0; l < QZ_CROSS; l++) 
  for (int m = 0; m < QZ_VAR; m++) 
  for (int n = 0; n < QZ_BVAR; n++) 
    symb_chs_init(&ctx->diff[i][c][j][k][l][m][n]);

/*
  for (int i = 0; i < QZ_SIZE_S; i++) 
  for (int j = 0; j < 8; j++) 
  for (int k = 0; k < QZ_BVAR_S; k++) 
  for (int l = 0; l < QZ_BVAR_S; l++) 
  for (int m = 0; m < QZ_BVAR_S; m++) 
    symb_chs_init(&ctx->mask[j][i][k][l][m]);
*/

  symb_chs_init(&ctx->methods);

  for (int c = 0; c < 3; c++) {
    symb_chs_init(&ctx->r_min[c]);
    symb_chs_init(&ctx->r_max[c]);
  }
}

// calculate color distance between two pixels
void static inline color_distance(pixel_t *p1, pixel_t *p2, int* dist) {
  dist[0] = COLOR_DIST(p1->d[0],p2->d[0]);
  dist[1] = COLOR_DIST(p1->d[1],p2->d[1]);
  dist[2] = COLOR_DIST(p1->d[2],p2->d[2]);
  return;
}

uint32_t static image_crc(image_t *i) {
  uint_fast32_t crc=0;
  crc32k_transform(crc,i->w & 255);
  crc32k_transform(crc,i->w / 256);
  crc32k_transform(crc,i->h & 255);
  crc32k_transform(crc,i->h / 256);
  for (int y = 0; y < i->h; y++) {
    for (int x = 0; x < i->w; x++) {
      pixel_t *p = image_pixel(i,x,y);
      for (int c = 0; c < 3; c++) {
        crc32k_transform(crc,p->d[c] & 255);
        crc32k_transform(crc,p->d[c] / 256);
      }
    }
  }
  return (~crc & 0xFFFFFFFF);
}

void output_interpolation_methods(coder_t *encode) {
        symb_put_int(&encode->coder,
                     &encode->ctx->methods,
                     encode->nb_interpolation_methods-1,
                     0,
                     MAX_INTERPOLATIONS_EVER,
                     NULL
                    );
        for(int i=0; i<encode->nb_interpolation_methods; i++) {
                symb_put_int(&encode->coder,
                     &encode->ctx->methods,
                     encode->interpolation_methods[i],
                     0,
                     MAX_INTERPOLATIONS_EVER,
                     NULL
                );
        }
}
void input_interpolation_methods(coder_t *decode) {
        decode->nb_interpolation_methods
        = symb_get_int(&decode->coder,
                       &decode->ctx->methods,
                       0,
                       MAX_INTERPOLATIONS_EVER
                      ) + 1;
        if (decode->nb_interpolation_methods > MAX_INTERPOLATIONS) {
                fprintf(stderr,"WARNING: using more interpolation methods (%i) than supported by this implementation (%i).\n",decode->nb_interpolation_methods, MAX_INTERPOLATIONS);
        }
        for(int i=0; i<decode->nb_interpolation_methods; i++) {
                decode->interpolation_methods[i] 
                = symb_get_int(&decode->coder,
                     &decode->ctx->methods,
                     0,
                     MAX_INTERPOLATIONS_EVER
                );
                if (decode->interpolation_methods[i] > MAX_INTERPOLATIONS-1) {
                        fprintf(stderr,"WARNING: using interpolation method %i, which is not supported by this implementation.\n",decode->interpolation_methods[i]);
                }
        }        
}
void output_interpolation_choice(coder_t *encode, int i, int c, int method, int prev) {
     int j=0;
     while( j<method ) {
//        symb_put_bit(&encode->coder,&encode->ctx->interpolationMethod[i][c][j][prev],1);
        symb_put_bit(&encode->coder,&encode->ctx->interpolationMethod[i][c][j][prev],1,&encode->ctx->bitsInterpol[c]);
        encode->ctx->symbInterpol[c]++;
        j++;
     }
//     if(j < encode->nb_interpolation_methods-1) symb_put_bit(&encode->coder,&encode->ctx->interpolationMethod[i][c][j][prev],0);
     if(j < encode->nb_interpolation_methods-1) {
       symb_put_bit(&encode->coder,&encode->ctx->interpolationMethod[i][c][j][prev],0,&encode->ctx->bitsInterpol[c]);
       encode->ctx->symbInterpol[c]++;
     }
}

int input_interpolation_choice(coder_t *decode, int i, int c, int prev) {
     int j=0;
     while ( j < decode->nb_interpolation_methods-1 
//             && symb_get_bit(&decode->coder,&decode->ctx->interpolationMethod[i][c][j][prev])) {
             && symb_get_bit(&decode->coder,&decode->ctx->interpolationMethod[i][c][j][prev])) {
        j++;
     }
     
     return j;
}

void output_mask(coder_t *encode,int oldmask,int newmask,int qsize,uint32_t *bvar) {
        assert(oldmask <= newmask);
//        int maskdiff = newmask-oldmask;
//        assert(maskdiff <= 7-oldmask);
        assert(qsize < QZ_SIZE_S);
/*        
        symb_put_int(&encode->coder,
                     &encode->ctx->mask[oldmask][qsize][qbvar(bvar[0],QZ_BVAR_S)][qbvar(bvar[1]+bvar[2],QZ_BVAR_S)][0*qbvar(bvar[2],QZ_BVAR_S)],
                     maskdiff,
                     0,
                     7-oldmask,
                     NULL
                    );
                    */
//    encode->ctx->symbSplit[3]++;
//    if (oldmask == newmask) {
//     symb_put_bit(&encode->coder,&encode->ctx->splitContA[qsize][oldmask][qbvar(bvar[0]+bvar[1]+bvar[2],QZ_BVAR_S)],0,&encode->ctx->bitsSplit[3]);
//    } else {
//     symb_put_bit(&encode->coder,&encode->ctx->splitContA[qsize][oldmask][qbvar(bvar[0]+bvar[1]+bvar[2],QZ_BVAR_S)],1,&encode->ctx->bitsSplit[3]);
     for (int c = 0; c < 3 ; c++) {
      if (!(oldmask & (1 << c))) {
       int others=0;
       if (c==0) others = oldmask/2;
       if (c==1) others = ((oldmask & 4)>>1) | (newmask & 1);
       if (c==2) others = newmask & 3;
       if (others< 0 || others > 3) fprintf(stderr,"oei %i\n",others);
       if ((oldmask & (1 << c)) != (newmask & (1 << c))) {
//        if (c==0 && ((oldmask&6) == (newmask&6))) {
//          // do nothing, it is implied that plane 0 had a split
//        } else {
          symb_put_bit(&encode->coder,&encode->ctx->splitCont[qsize][c][others][qbvar(bvar[c],QZ_BVAR_S)],1,&encode->ctx->bitsSplit[c]);
          encode->ctx->symbSplit[c]++;
//        }
       } else {
          symb_put_bit(&encode->coder,&encode->ctx->splitCont[qsize][c][others][qbvar(bvar[c],QZ_BVAR_S)],0,&encode->ctx->bitsSplit[c]);
          encode->ctx->symbSplit[c]++;
          encode->splits[c]++;
       }
     }
    }
//   }
}

int input_mask(coder_t *decode,int mask,int qsize,uint32_t *bvar) {
     int oldmask = mask;
     int newmask = oldmask;
     for (int c = 0; c < 3 ; c++) {
      if (!(oldmask & (1 << c))) {
       int others=0;
       if (c==0) others = oldmask/2;
       if (c==1) others = ((oldmask & 4)>>1) | (newmask & 1);
       if (c==2) others = newmask & 3;
       if (symb_get_bit(&decode->coder, &decode->ctx->splitCont[qsize][c][others][qbvar(bvar[c],QZ_BVAR_S)] )) {
                newmask |= 1<<c;
       }
     }
    }
    return newmask;
}

void output_min_max(coder_t *encode,int min,int max,int oldmin,int oldmax, int c) {
        assert(min <= max);
        symb_put_int(&encode->coder,
                     &encode->ctx->r_min[c],
                     min,
                     oldmin,
                     oldmax,
                     NULL
                    );
        symb_put_int(&encode->coder,
                     &encode->ctx->r_max[c],
                     (oldmax-max),
                     0,
                     (oldmax-min),
                     NULL
                    );
//        fprintf(stderr,"min: %i  (range: %i to %i)\n",min,oldmin,oldmax);                    
  //      fprintf(stderr,"max: %i  (range: %i to %i)\n",oldmax-max,0,oldmax-min);                    
    //    fprintf(stderr,"OUTPUT %i,%i,%i,%i,%i \n",min,max,oldmin,oldmax,c);                    
}
void input_min_max(coder_t *decode,int *min,int *max,int oldmin,int oldmax, int c) {
        *min = symb_get_int(&decode->coder,
                     &decode->ctx->r_min[c],
                     oldmin,
                     oldmax
                    );
        *max = oldmax - symb_get_int(&decode->coder,
                     &decode->ctx->r_max[c],
                     0,
                     (oldmax - *min)
                    );
//        fprintf(stderr,"min: %i   max: %i\n",*min,*max);                    
  //      fprintf(stderr,"min: %i  (range: %i to %i)\n",*min,oldmin,oldmax);                    
    //    fprintf(stderr,"max: %i  (range: %i to %i)\n",oldmax-*max,0,oldmax-*min);                    
      //  fprintf(stderr,"INPUT %i,%i,%i,%i,%i \n",*min,*max,oldmin,oldmax,c);                    
}



// output the color values (from planes not given by mask) of a pixel to the range encoder
void static output_pixel(coder_t *encode, int x, int y, const pixel_t *guess, uint32_t size, int mask, 
                                const pixel_t *topleft, const pixel_t *topright, const pixel_t *left,const pixel_t *right, uint32_t bvar[3]) {
//  const pixel_t def = {.d= {STRETCH(128),STRETCH(256),STRETCH(256)}};

//  fprintf(stderr,"pixel (%i,%i)...\n",x,y);

  const pixel_t def = {.d= {128,256,256}};
  pixel_t *p = image_pixel(encode->source,x,y);
  pixel_t *r = image_pixel(encode->rec,x,y);
//        assert(r->d[2] != path[2]);
  //  pixel_t *r=image_pixel(encode->rec,x,y);
  if (!guess) guess = &def;
  if (!topleft) topleft = &def;
  if (!topright) topright = &def;
  if (!left) left = &def;
  if (!right) right = &def;
    for (int c = 0; c < 3; c++) {
      if (!(mask & (1 << c))) {
        encode->outputted_pixels[c]++;
#if (ACCURATE_YIQ_RANGES == 1)
        int minval=0;
        int maxval=0;

        switch (c) {
          case 0:
            get_range_y(&minval,&maxval,encode->range_min,encode->range_max);
            break;
          case 1:
            get_range_i(r->d[0],&minval,&maxval,encode->range_min,encode->range_max);
            break;
          case 2:
            get_range_q(r->d[0],r->d[1],&minval,&maxval,encode->range_min,encode->range_max);
            break;
        }
        if (p->d[c] < minval) p->d[c] = minval;
        if (p->d[c] > maxval) p->d[c] = maxval;
#else
//        int minval=STRETCH(0),maxval=STRETCH(c==0 ? 255 : 510);
        int minval = STRETCH(encode->range_min[c]);
        int maxval = STRETCH(encode->range_max[c]);
#endif
        int qz = encode->qz[c];
        int dist = -rdiv(guess->d[c] - (p->d[c]-(qz/2)),qz);
        int min = -rdiv(guess->d[c] - (minval-(qz/2)),qz);
        int max = -rdiv(guess->d[c] - (maxval-(qz/2)),qz);

//        fprintf(stderr,"pixel %i %c(%i,%i)\n",c,planesYIQ[c],x,y);
        
        int lrdiff = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]),(c==0 ? 255 : 510),0,QZ_LRDIFF);
        int grad   = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]), (c==0 ? 510 : 1020),0,QZ_GRAD);
        int cross  = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(right->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(topright->d[c]), (c==0 ? 510 : 1020),0,QZ_CROSS);
        int variance  = quantize_log_uint32(
                        (c==0 ? 4000 : 1000) *
                        (4*UNSTRETCH(topleft->d[c])*UNSTRETCH(topleft->d[c])
                        +4*UNSTRETCH(topright->d[c])*UNSTRETCH(topright->d[c])
                        +4*UNSTRETCH(left->d[c])*UNSTRETCH(left->d[c])
                        +4*UNSTRETCH(right->d[c])*UNSTRETCH(right->d[c])
                -
                        UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c])
                      * UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c])
                        ),
                        QZ_VAR);

        int qsize  = quantize_log_uint32(size,QZ_SIZE+1)-1;
        
        
        assert(lrdiff>=0);
        assert(grad>=0);
        assert(cross>=0);
        assert(lrdiff<QZ_LRDIFF);
        assert(grad<QZ_GRAD);
        assert(cross<QZ_CROSS);
        assert(qsize>=0);
        assert(variance>=0);

//        fprintf(stderr,"context diff[%i][%i][%i][%i][%i] for value %i in range [%i,%i]\n",qsize,c,lrdiff,grad,cross,dist,min,max);
        
        symb_put_int(&encode->coder,
                     &encode->ctx->diff[qsize][c][lrdiff][grad][cross][variance][qbvar(bvar[c],QZ_BVAR)],
                     dist,min,max,
                     &encode->ctx->bitsPixeldata[c]
                    );
        encode->ctx->symbPixeldata[c]++;

        r->d[c] = dist*qz + guess->d[c];


/*  //uncomment to see trajectory
        r->d[c] = path[c];
        pixelcounter++;
        if (pixelcounter%90 == 0) {path[0]++;}
  */      
        
//        if (pixelcounter%50 == 0) {path[1]++;}
        
//        assert(encode->maxD[c]==0 ? (COLOR_DIST(r->d[c],p->d[c])==0) : 1);
      }
    }
//   fprintf(stderr,"pixel (%i,%i) done\n",x,y);
//    pixel_set(r,p,mask);
    //    fprintf(stderr,"output diff: %i %i %i\n",p->d[0]-guess->d[0],p->d[1]-guess->d[1],p->d[2]-guess->d[2]);
}

// input the color values (from planes not given by mask) of a pixel from the range encoder
void static input_pixel(coder_t *decode, int x, int y, const pixel_t *guess, uint32_t size, int mask,
                                const pixel_t *topleft, const pixel_t *topright, const pixel_t *left,const pixel_t *right, uint32_t bvar[3]) {
  const pixel_t def = {.d= {128,256,256}};
  //pixel_t *p = image_pixel(decode->source,x,y);
  pixel_t *r = image_pixel(decode->rec,x,y);
  if (!guess) guess = &def;
  if (!topleft) topleft = &def;
  if (!topright) topright = &def;
  if (!left) left = &def;
  if (!right) right = &def;
    for (int c = 0; c < 3; c++) {
      if (!(mask & (1 << c))) {
#if (ACCURATE_YIQ_RANGES == 1)
        int minval=0, maxval=0;
        switch (c) {
          case 0:
            get_range_y(&minval,&maxval,decode->range_min,decode->range_max);
            break;
          case 1:
            get_range_i(r->d[0],&minval,&maxval,decode->range_min,decode->range_max);
            break;
          case 2:
            get_range_q(r->d[0],r->d[1],&minval,&maxval,decode->range_min,decode->range_max);
            break;
        }
#else
//        int minval=STRETCH(0),maxval=STRETCH(c==0 ? 255 : 510);
        int minval = STRETCH(decode->range_min[c]);
        int maxval = STRETCH(decode->range_max[c]);
#endif
        int qz = decode->qz[c];
        int min = -rdiv(guess->d[c] - (minval-(qz/2)),qz);
        int max = -rdiv(guess->d[c] - (maxval-(qz/2)),qz);

        int lrdiff = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]),(c==0 ? 255 : 510),0,QZ_LRDIFF);
        int grad   = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(right->d[c]), (c==0 ? 510 : 1020),0,QZ_GRAD);
        int cross  = quantize_log(UNSTRETCH(topleft->d[c])+UNSTRETCH(right->d[c])-UNSTRETCH(left->d[c])-UNSTRETCH(topright->d[c]), (c==0 ? 510 : 1020),0,QZ_CROSS);
        int variance  = quantize_log_uint32(
                        (c==0 ? 4000 : 1000) *
                        (4*UNSTRETCH(topleft->d[c])*UNSTRETCH(topleft->d[c])
                        +4*UNSTRETCH(topright->d[c])*UNSTRETCH(topright->d[c])
                        +4*UNSTRETCH(left->d[c])*UNSTRETCH(left->d[c])
                        +4*UNSTRETCH(right->d[c])*UNSTRETCH(right->d[c])
                -
                        UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c])
                      * UNSTRETCH(topleft->d[c])+UNSTRETCH(topright->d[c])+UNSTRETCH(left->d[c])+UNSTRETCH(right->d[c])
                        ),
                        QZ_VAR);

        int qsize  = quantize_log_uint32(size,QZ_SIZE+1)-1;

        assert(qsize>=0);

        int dist = symb_get_int(&decode->coder,&decode->ctx->diff[qsize][c][lrdiff][grad][cross][variance][qbvar(bvar[c],QZ_BVAR)],min,max);
        r->d[c] = dist*qz + guess->d[c];
      }
    }
}


// calculate color distance between an area and its interpolation
void static area_distance(coder_t *encode, int x1, int x2, int y1, int y2, int *maxd, double *avgd, uint64_t *ndiff) {
  uint64_t fd[3] = {0,0,0};
  double qsumD = 0.0;
  int maxD = 0, maxDi = 0, maxDq = 0;
  for (int y = y1; y <= y2; y++) {
    for (int x = x1; x <= x2; x++) {
      pixel_t *guess = image_pixel(encode->rec,x,y);
      pixel_t *real = image_pixel(encode->source,x,y);
      int distance[3];
      color_distance(real,guess,distance);
      fd[0] += distance[0]; //(real->d[0]-guess.d[0])*(real->d[0]-guess.d[0]);
      fd[1] += distance[1]; //(real->d[1]-guess.d[1])*(real->d[1]-guess.d[1]);
      fd[2] += distance[2]; //(real->d[2]-guess.d[2])*(real->d[2]-guess.d[2]);
      if (distance[0] > maxD) maxD = distance[0];
      if (distance[1] > maxDi) maxDi = distance[1];
      if (distance[2] > maxDq) maxDq = distance[2];
      qsumD += distance[0] * distance[0];
    }
  }
  uint32_t pixels = ((y2 - y1 + 1) * (x2 - x1 + 1));
  maxd[0] = maxD;
  maxd[1] = maxDi;
  maxd[2] = maxDq;
  *avgd = sqrt(qsumD / pixels);
  ndiff[0] = fd[0];
  ndiff[1] = fd[1];
  ndiff[2] = fd[2];
}

// recursion types (not implemented yet)
typedef enum {
  RECURSE_LONG, RECURSE_ASSYM_LEFT, RECURSE_ASSYM_RIGHT, RECURSE_CENTER,
} recurse_t;


// order of guess: left, right, top, bottom
void static border_guess(int guess_method, coder_t *coder, int px1, int px2, int py1, int py2, int xm, int ym, pixel_t *guess, int mask) {
//  fprintf(stderr,"border_guess(%i,%i - %i,%i - %i,%i
  switch(guess_method) {
    case 0:
          pixel_linear(&guess[0],image_pixel(coder->rec,px1,py1),ym - py1,image_pixel(coder->rec,px1,py2),py2 - ym,mask);
          pixel_linear(&guess[1],image_pixel(coder->rec,px2,py1),ym - py1,image_pixel(coder->rec,px2,py2),py2 - ym,mask);
          pixel_linear(&guess[2],image_pixel(coder->rec,px1,py1),xm - px1,image_pixel(coder->rec,px2,py1),px2 - xm,mask);
          pixel_linear(&guess[3],image_pixel(coder->rec,px1,py2),xm - px1,image_pixel(coder->rec,px2,py2),px2 - xm,mask);
          return;

    case 1:
          gradient(&guess[0],coder->rec,px1,ym,px1,px2,py1,py2,mask);
          gradient(&guess[1],coder->rec,px2,ym,px1,px2,py1,py2,mask);
          gradient(&guess[2],coder->rec,xm,py1,px1,px2,py1,py2,mask);
          gradient(&guess[3],coder->rec,xm,py2,px1,px2,py1,py2,mask);
          return;


    case 2:
          guess[0]=purple;
          guess[1]=purple;
          guess[2]=purple;
          guess[3]=purple;
          return;
  }
}

int static no_interpolation_choice(coder_t *encode, int x1, int x2, int y1, int y2, int mask) {
/*       
        if (pixel_equal(image_pixel(encode->rec,x1,y1),image_pixel(encode->rec,x1,y2),mask)
         && pixel_equal(image_pixel(encode->rec,x1,y1),image_pixel(encode->rec,x2,y1),mask)
         && pixel_equal(image_pixel(encode->rec,x1,y1),image_pixel(encode->rec,x2,y2),mask)
           ) return 1;
*/
        return 0;
}


void static compute_bvar(coder_t *encode, int x1, int y1, int x2, int y2, int mask, uint32_t *bvar, int hastop, int hasright, int hasbottom, int hasleft) {
  return;
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  int val;
  for(int c=0;c<3;c++) {
   if(!(mask & 1<<c)) {
    uint64_t sum=0;
    uint64_t sumq=0;
    uint64_t count=0;
    if (hasleft) {
     for(int y=y1; y<y2; y++) {
         val=UNSTRETCH(image_pixel(encode->rec,px1,y)->d[c]);
         sum += val;
         sumq += val*val;
         count++;
     }
    }
    if (hasright) {
     for(int y=y1; y<y2; y++) {
         val=UNSTRETCH(image_pixel(encode->rec,px2,y)->d[c]);
         sum += val;
         sumq += val*val;
         count++;
     }
    }
    if (hastop) {
     for(int x=x1; x<x2; x++) {
         val=UNSTRETCH(image_pixel(encode->rec,x,py1)->d[c]);
         sum += val;
         sumq += val*val;
         count++;
     }
    }
    if (hasbottom) {
     for(int x=x1; x<x2; x++) {
         val=UNSTRETCH(image_pixel(encode->rec,x,py2)->d[c]);
         sum += val;
         sumq += val*val;
         count++;
     }
    }
    int mult=8;         // weight of known corners (4 = same as other known border pixels, 8 = counted twice, etc.)

    if (count<10) mult=4;
//    if (count>30) mult=4;

//    fprintf(stderr,"\ncount:%i\n",count);
    val=UNSTRETCH(image_pixel(encode->rec,px1,py1)->d[c]);
    sum += mult*val/4;
    sumq += mult*val*val/4;
    val=UNSTRETCH(image_pixel(encode->rec,px1,py2)->d[c]);
    sum += mult*val/4;
    sumq += mult*val*val/4;
    val=UNSTRETCH(image_pixel(encode->rec,px2,py1)->d[c]);
    sum += mult*val/4;
    sumq += mult*val*val/4;
    val=UNSTRETCH(image_pixel(encode->rec,px2,y2)->d[c]);
    sum += mult*val/4;
    sumq += mult*val*val/4;
    count += mult;
    if (count>0) {
      int scale = 16000;
      if (c>0) scale = 64000;
      uint32_t var=(sumq*scale/count)-(sum*sum*scale/count/count);
      if (sumq < sum*sum/count) fprintf(stderr,"0AAAOOOOOOOOOOOOOOOOO (count %llu) %llu < %llu ^2 /count\n",(unsigned long long)count,(unsigned long long)sumq,(unsigned long long)sum);
      sum *= sum*scale;
      sum /= count;
      sumq *= scale;
      if (sumq < sum) fprintf(stderr,"1AAAOOOOOOOOOOOOOOOOO (count %llu) %llu < %llu\n",(unsigned long long)count,(unsigned long long)sumq,(unsigned long long)sum);
      sum /= count;
      sumq /= count;
      uint32_t var3= sumq - sum;
      if (var != var3) fprintf(stderr,"OOOOOOOOOOOOOOOOO %i != %i\n",var,var3);
      if (sumq < sum) fprintf(stderr,"2AAAOOOOOOOOOOOOOOOOO (count %llu) %llu < %llu\n",(unsigned long long)count,(unsigned long long)sumq,(unsigned long long)sum);
//      long long int var2=(sumq*scale/count)-(sum*sum*scale/count/count);
//      if (var2<0) fprintf(stderr,"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE %lli = %llu - %llu\n",var2,(sumq*scale/count),(sum*sum*scale/count/count));
      bvar[c] = var;
      //quantize_log_uint32(2*var+10,quant+4)-4;
//      fprintf(stderr,"[%i] count:%lli   qbvar:%i   bvar:%i\n",c,count,bvar[c],var);
//      fprintf(stderr,"v:%i\t",bvar[c]);
      assert(bvar[c]>=0);
      if (bvar[c] < 0) bvar[c]=0;
    } else {
      bvar[c] = 0;
    }
    stats_bvar[c][qbvar(bvar[c],QZ_BVAR)]++;
   }
  }
}

int static nn(int x, int y, int px1, int py1, int px2, int py2) {
        if (x == px1 && y == py1) return 0;
        if (x == px1 && y == py2) return 0;
        if (x == px2 && y == py1) return 0;
        if (x == px2 && y == py2) return 0;
        return 1;
}

// recursive function for encoding
void static encode_recurse(coder_t *encode, int x1, int x2, int y1, int y2, int mask, int previous_method, int hastop, int hasright, int hasbottom, int hasleft, int shape, int dir) {
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  uint32_t p = w * h;
//  fprintf(stderr,"w:%i  h:%i   area:%i\n",w,h,p);
  int size = ilog2(p);
  int qsize = quantize_log_uint32(p,QZ_SIZE_S+1)-1;
//  int qsize = (size<=QZ_SIZE_S/2 ? size : QZ_SIZE_S/2 + quantize_log_uint32(p,QZ_SIZE_S/2+1)-1);
//  fprintf(stderr,"entering (%i,%i)-(%i,%i)\n",x1,y1,x2,y2);
  
  
  if (x2 < x1 || y2 < y1) return;

//  if (px2 - px1 <= 1 && py2 - py1 <= 1) {
//    interpolation(3,encode->rec,x1,x2,y1,y2,mask);
//    return;
//  }

  if (px2 - px1 <= 1 && py2 - py1 <= 1) return;
  
//  if (x2 - x1 == 0 && y2 - y1 == 0) {
//  }

/*
  fprintf(stderr,"area (%i,%i)-(%i,%i) ",x1,y1,x2,y2);
  fprintf(stderr,"dir:%i shape:%i ",dir,shape);
  if (hastop) fprintf(stderr,"top ");
  if (hasright) fprintf(stderr,"right ");
  if (hasbottom) fprintf(stderr,"bottom ");
  if (hasleft) fprintf(stderr,"left ");
  fprintf(stderr,"\n");
*/

  int method;
  int current_method = previous_method;

  uint32_t bvar[3] = {};
  compute_bvar(encode,x1,y1,x2,y2,mask,&bvar[0],hastop,hasright,hasbottom,hasleft);


  if (p <= encode->minsplit) {
    double avgd;
    int maxd[3];
    uint64_t rd[3];
    int best_maxd[3]={};
    uint64_t best_rd[3]={};;
    uint64_t second_best_rd[3]={-1,-1,-1};;
    int maxd_i[3]={-1,-1,-1};
    int rd_i[3]={-1,-1,-1};
    int choice[3]={1,1,1};
    
    for(int i=0; i<encode->nb_interpolation_methods; i++) {
      interpolation(encode->interpolation_methods[i],encode->rec,x1,x2,y1,y2,mask,hastop,hasright,hasbottom,hasleft);
      area_distance(encode,x1,x2,y1,y2,maxd,&avgd,rd);
      for (int c = 0; c < 3; c++) {
        if (!(mask & (1 << c)) && choice[c]) {
//        fprintf(stderr,"Method: %i channel %i, maxd=%i\n",encode->interpolation_methods[i],c,maxd[c]);
          if (no_interpolation_choice(encode,x1,x2,y1,y2,~(1<<c))) {
            best_maxd[c]=maxd[c]; maxd_i[c]=i;
            best_rd[c]=rd[c]; rd_i[c]=i; 
            choice[c]=0;
          } else {
            if (maxd_i[c] == -1 || best_maxd[c] > maxd[c]) { best_maxd[c]=maxd[c]; maxd_i[c]=i; }
            if (rd_i[c] == -1 || best_rd[c] > rd[c]) { 
              second_best_rd[c] = best_rd[c]; best_rd[c]=rd[c]; rd_i[c]=i; 
            } else {
              if (second_best_rd[c] == -1 || second_best_rd[c] > rd[c]) second_best_rd[c] = rd[c];
            }
          }
        }
      }
    }
    for (int c = 0; c < 3; c++) {maxd[c] = best_maxd[c]; rd[c] = best_rd[c];}



    int newmask = mask;
    double factor = encode->factor;

    // restrict factor in smaller areas:  
    if (factor > 3 * p) factor = 3 * p;

    // uncomment to have average error factor
    // factor *= p;

    // stop criterion
    for (int c = 0; c < 3 ; c++) {
      if (!(mask & (1 << c))) {
        method=maxd_i[c];
//        fprintf(stderr,"[%i] newmask %i   others %i\n",c,newmask,others);
        if (   (maxd[c] <= encode->maxD[c] )   || rd[c] <= encode->maxD[c] * factor ) {
          if (maxd[c] > encode->maxD[c]) method=rd_i[c];
          if (maxd[c] > maxdist[c]) maxdist[c] = maxd[c];
          interpolation(encode->interpolation_methods[method],encode->rec,x1,x2,y1,y2,~(1<<c),hastop,hasright,hasbottom,hasleft);
          if (choice[c]) {
            encode->method[encode->interpolation_methods[method]]++;
            encode->method_gain[encode->interpolation_methods[method]] += second_best_rd[c]-rd[c];
          } else {
            encode->method[MAX_INTERPOLATIONS]++;
          }
          newmask |= (1 << c);
        }
      }
    }

    output_mask(encode,mask,newmask,qsize,bvar);

     for (int c = 0; c < 3 ; c++) {
      if (!(mask & (1 << c))) {
        method=maxd_i[c];
        if ((mask & (1 << c)) != (newmask & (1 << c))) {
          output_interpolation_choice(encode,size,c,method,current_method);
        }
      }
     }
    
    
    if (newmask == 7) return;

    mask = newmask;
  //  fprintf(stderr,"real dist: %lu %lu %lu\n",(unsigned long)(rd[0]),(unsigned long)(rd[1]),(unsigned long)(rd[2]));
  }


//  fprintf(stderr,"splitting (%i,%i)-(%i,%i)\n",x1,y1,x2,y2);

  pixel_t *orig[4] = {
        image_pixel(encode->rec,px1,py1),               //left top
        image_pixel(encode->rec,px2,py1),               //right top
        image_pixel(encode->rec,px1,py2),               //left bottom
        image_pixel(encode->rec,px2,py2)                //right bottom
        };

  pixel_t guess[4];
// dir=1;
  int xm = (px1 + px2) / 2 ;
  int ym = (py1 + py2) / 2 ;
  border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,mask);


 
//  fprintf(stderr,"recursive calls\n");


  if (px2 - px1 > py2 - py1 && shape!=1) {
//  if (x2 - x1 > y2 - y1) {
//  if (!shape) {
    if (!hastop && nn(xm,y1,px1,py1,px2,py2)) { encode->pTop++; output_pixel(encode,xm,y1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar); }
    if (!hasbottom && nn(xm,y2,px1,py1,px2,py2) && (y1 != y2 || !(!hastop && nn(xm,y1,px1,py1,px2,py2)))) { encode->pBottom++; output_pixel(encode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar);}
    if (shape==2) {
    if (dir == 0) {
      encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   0,0); // left side
      encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,1,0); // right side
    } else if (dir == 1) {
      encode_recurse(encode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, 1,0); // right side
      encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,0,1); // left side
    } else if (dir == 2) {
      encode_recurse(encode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, 0,1); // right side
      encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,1,1); // left side
    } else if (dir == 3) {
      encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   1,1); // left side
      encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,0,0); // right side
    }
    } else {
    if (dir == 0) {
      encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   0,0); // left side
      encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,0,1); // right side
    } else {
      encode_recurse(encode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, 0,1); // right side
      encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,0,0); // left side
    }
    }
  } else {
    if (!hasleft && nn(x1,ym,px1,py1,px2,py2)) { encode->pLeft++; output_pixel(encode,x1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar); }
    if (!hasright && nn(x2,ym,px1,py1,px2,py2) && (x1!=x2 || !(!hasleft && nn(x1,ym,px1,py1,px2,py2)))) { encode->pRight++; output_pixel(encode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar); }
    if (shape==1) {
    if (dir == 0) {
      encode_recurse(encode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,     2,0); // top side
      encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,2,1); // bottom side
    } else {
      encode_recurse(encode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft,  2,2); // bottom side
      encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   2,3); // top side
    }
    } else {
    if (dir == 0) {
      encode_recurse(encode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,     1,0); // top side
      encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,0,0); // bottom side
    } else {
      encode_recurse(encode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft,  0,0); // bottom side
      encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   1,1); // top side
    }
    }
  }
//  fprintf(stderr,"done area (%i,%i)-(%i,%i)\n",x1,y1,x2,y2);
}

// recursive function for decoding
void static decode_recurse(coder_t *decode, int x1, int x2, int y1, int y2, int mask, int previous_method, int hastop, int hasright, int hasbottom, int hasleft, int shape, int dir) {
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  uint32_t p = w * h;
  int size = ilog2(p);
  int qsize = quantize_log_uint32(p,QZ_SIZE_S+1)-1;
  int method;
  int current_method = previous_method;

//  if (px2 - px1 <= 1 && py2 - py1 <= 1) {
  if (x2 - x1 == 0 && y2 - y1 == 0) {
//    interpolate(decode->rec,x1,x2,y1,y2,mask);
    return;
  }

  uint32_t bvar[3] = {};
  compute_bvar(decode,x1,y1,x2,y2,mask,&bvar[0],hastop,hasright,hasbottom,hasleft);

  int oldmask = mask;
  if (p <= decode->minsplit) {
  
     mask = input_mask(decode,oldmask,qsize,bvar);
//    mask=7;
    for (int c = 0; c < 3; c++) {
      if (!(oldmask & (1 << c))) {
        if (no_interpolation_choice(decode,x1,x2,y1,y2,~(1<<c))) {
           int i=0;
           interpolation(decode->interpolation_methods[i],decode->rec,x1,x2,y1,y2,~(1<<c),hastop,hasright,hasbottom,hasleft);
        } else {
          if (mask & (1 << c)) {
            method = input_interpolation_choice(decode,size,c,current_method);
            interpolation(decode->interpolation_methods[method],decode->rec,x1,x2,y1,y2,~(1<<c),hastop,hasright,hasbottom,hasleft);
            decode->method[decode->interpolation_methods[method]]++;
          }
        }
      }
    }
/*
    if (oldmask != mask) {
          interpolate(decode->rec,x1,x2,y1,y2,oldmask | ~mask);
    }
*/
    if (mask == 7) return;
  }

  /*  recurse_t rec=RECURSE_CENTER;
   if (w > h) {
   if (h >= MIN_LONG_WIDTH && p >= MIN_LONG_PIXELS) {
   if (jif_get_bit(decode,&decode->ctx->splitShort[size][c])) {
   if (w >= MIN_ASSYM_WIDTH && p >= MIN_ASSYM_PIXELS) {
   if (jif_get_bit(decode,&decode->ctx->splitAssym[size][c])) {
   if (jif_get_bit(decode,&decode->ctx->splitSided[size][c])) {

   } else {
   }
   } else {
   rec=RECURSE_CENTER;
   }
   }
   } else {
   rec=RECURSE_LONG;
   }
   }
   }*/

  pixel_t *orig[4] = {
        image_pixel(decode->rec,px1,py1),               //left top
        image_pixel(decode->rec,px2,py1),               //right top
        image_pixel(decode->rec,px1,py2),               //left bottom
        image_pixel(decode->rec,px2,py2)                //right bottom
        };
  pixel_t guess[4];
  int xm = (x1 + x2) / 2;
  int ym = (y1 + y2) / 2;
  border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,mask);

//  dir=0;

/*
  if (px2 - px1 > py2 - py1) {
    int xm = (px1 + px2) / 2;
    if (y1==0) input_pixel(decode,xm,py1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar);
    input_pixel(decode,xm,py2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar);
    decode_recurse(decode,x1,xm,y1,y2,mask,current_method); // left side
    decode_recurse(decode,xm + 1,x2,y1,y2,mask,current_method); // right side
  } else {
    int ym = (py1 + py2) / 2;
    if (x1==0) input_pixel(decode,px1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar);
    input_pixel(decode,px2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar);
    decode_recurse(decode,x1,x2,y1,ym,mask,current_method); // top side
    decode_recurse(decode,x1,x2,ym + 1,y2,mask,current_method); // bottom side
  }*/
//  if (!shape) {
//  if (px2 - px1 > py2 - py1) {
  if (x2 - x1 > y2 - y1) {
    if (!hastop) { decode->pTop++; input_pixel(decode,xm,y1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar); }
    if (!hasbottom) { decode->pBottom++; input_pixel(decode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar);}
    if (dir == 0) {
      decode_recurse(decode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,!shape,0); // left side
      decode_recurse(decode,xm + 1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,!shape,1); // right side
    } else {
      decode_recurse(decode,xm + 1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,!shape,1); // right side
      decode_recurse(decode,x1,xm,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,!shape,0); // left side
    }
  } else {
    if (!hasleft) { decode->pLeft++; input_pixel(decode,x1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar); }
    if (!hasright) { decode->pRight++; input_pixel(decode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar); }
    if (dir == 0) {
      decode_recurse(decode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,!shape,0); // top side
      decode_recurse(decode,x1,x2,ym + 1,y2,mask,current_method,1,hasright,hasbottom,hasleft,!shape,1); // bottom side
    } else {
      decode_recurse(decode,x1,x2,ym + 1,y2,mask,current_method,0,hasright,hasbottom,hasleft,!shape,1); // bottom side
      decode_recurse(decode,x1,x2,y1,ym,mask,current_method,hastop,hasright,1,hasleft,!shape,0); // top side
    }
  }
}


// decode an image
void static decode(image_t *img, char *in) {
  char c[4];
  FILE *f = fopen(in,"rb");
  int ret = fread(c,4,1,f);
  if (ret < 1 || memcmp(c,magic,4) != 0) {
    fprintf(stderr,"No JiF file\n");
    return;
  }
  int w = 0;
  int h = 0;
  w |= fgetc(f);
  w |= (fgetc(f) << 8);
  h |= fgetc(f);
  h |= (fgetc(f) << 8);
  int qzy = 0, qzi = 0, qzq = 0;
  qzy |= fgetc(f);
  qzy |= (fgetc(f) << 8);
  qzi |= fgetc(f);
  qzi |= (fgetc(f) << 8);
  qzq |= fgetc(f);
  qzq |= (fgetc(f) << 8);
  int cutoff = 0;
  cutoff |= fgetc(f);
  cutoff |= (fgetc(f) << 8) ;
  image_init(img,w,h,&purple);

  int min[3] = {0,0,0};
  int max[3] = {255,510,510};
/*  for(int c=0; c<3; c++) {
        min[c]=fgetc(f);
        max[c]=fgetc(f);
  }
  int msb = fgetc(f);
  min[1] |= (msb & 1<<0) << 8;
  max[1] |= (msb & 1<<1) << 7;
  min[2] |= (msb & 1<<2) << 6;
  max[2] |= (msb & 1<<3) << 5;
  */


  ctx_t *ctx=malloc(sizeof(ctx_t));
  coder_t dec = {.source = NULL,.rec = img, .qz={qzy,qzi,qzq}, .method={}, .method_gain={}, .ctx=ctx, .range_min={min[0],min[1],min[2]}, .range_max={max[0],max[1],max[2]}};
 
  fprintf(stderr,"Decoding %ix%i image with quantization (%f,%f,%f)\n",w,h,(double)qzy/COLOR_STRETCH,(double)qzi/COLOR_STRETCH,(double)qzq/COLOR_STRETCH);
  
  fprintf(stderr,"Using RAC range [%g-%g] (cutoff %i)\n",cutoff/4096.0,(4096-cutoff)/4096.0,cutoff);
  symb_init_read(&dec.coder,f,cutoff);
  ctx_init(dec.ctx,cutoff);
//  extra_ctx_init(&dec,cutoff);
  input_interpolation_methods(&dec);
  int skips=0;
  int ms=w*h;
  while (symb_get_simple_bit(&dec.coder)) { skips++; ms /= 2; }
  dec.minsplit=ms;
  fprintf(stderr,"Areas over %i pixels are always splitted\n",ms);
  image_set_alpha(img,0);
  int fail = 0;
  for (int i = 0; i < 16; i++) {
    int x = symb_get_simple_bit(&dec.coder);
    fail |= (x != ((imagic>>i) & 1));
  }
  if (fail) fprintf(stderr,"Header corrupt?\n");

  for(int c=0; c<3; c++) {
        input_min_max(&dec,&(dec.range_min[c]),&(dec.range_max[c]),0,(c==0?255:510),c);
        fprintf(stderr,"Channel %c: range=[%i,%i]\n",planesYIQ[c],dec.range_min[c],dec.range_max[c]);
  }


  uint32_t p= img->w * img->h;
  uint32_t bvar[3]={};
  input_pixel(&dec,0,0,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  input_pixel(&dec,img->w - 1,0,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  input_pixel(&dec,img->w - 1,img->h - 1,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  input_pixel(&dec,0,img->h - 1,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  decode_recurse(&dec,0,img->w - 1,0,img->h - 1,0,0,0,0,0,0,0,0);
  fail = 0;
  // 16 magic marker bits
  for (int i = 0; i < 16; i++) {
    int x = symb_get_simple_bit(&dec.coder);
    fail |= (x != ((imagic>>i) & 1));
  }
  // 32 crc bits
  uint32_t crc=image_crc(img);
  for (int i = 0; i < 32; i++) {
    int x = symb_get_simple_bit(&dec.coder);
    fail |= (x != ((crc>>i) & 1))*2;
  }
  if (fail & 1) {
    fprintf(stderr,"Image file is invalid\n");
  } else {
    if (fail & 2) {
      fprintf(stderr,"Checksum mismatch\n");
    } else {
      fprintf(stderr,"Checksum matches: 0x%08x\n",crc);
    }
  }
  yiq2rgb(img);
  fclose(f);
}

void parse_interpolation_methods(coder_t *enc, char *methods) {
        char * allowed="0123456789";
        int nb_methods = 0;
        while(methods[0]) {
            size_t length = strspn(methods,allowed);
            if (length==0) {methods++; continue;}
            char * end=NULL;
            unsigned long num = strtoul(methods,&end,10);
            if (end==methods) {methods++; continue; }   //TODO: give warning
            methods += length;
            if (num>MAX_INTERPOLATIONS-1) fprintf(stderr,"ERROR: I don't have an implementation of method %lu\n",num);
            enc->interpolation_methods[nb_methods]=num;
//            fprintf(stderr,"Method %i : %i\n",nb_methods,num);
            nb_methods++;
        }
        enc->nb_interpolation_methods=nb_methods;
        fprintf(stderr,"Nb of methods: %i \n",enc->nb_interpolation_methods);
}

// encode an image
void static encode(image_t *img, double epsilon, double epsiloni, double epsilonq, double factor, double qzy, double qzi, double qzq, char *out, char *methods, int skipsplit, int cutoff) {
  image_t rec = {};
  image_init(&rec,img->w,img->h,&purple);
  if (qzy < 1.0) qzy = 1.0;
  if (qzi < 1.0) qzi = 1.0;
  if (qzq < 1.0) qzq = 1.0;
  if (qzy > 255.0) qzy = 255.0;
  if (qzi > 255.0) qzi = 255.0;
  if (qzq > 255.0) qzq = 255.0;
  int ms=(img->w*img->h) >> skipsplit;
  ctx_t *ctx=malloc(sizeof(ctx_t));
  coder_t enc = {.source = img,.rec = &rec,.maxD = {epsilon,epsiloni,epsilonq},.factor = factor,
                 .outputted_pixels = {0,0,0}, .qz={(int)(qzy*COLOR_STRETCH+0.5), (int)(qzi*COLOR_STRETCH+0.5), (int)(qzq*COLOR_STRETCH+0.5)}, .pLeft=0, .pRight=0, .pTop=0, .pBottom=0,
                 .method={}, .method_gain={}, .minsplit=ms, .ctx=ctx, .range_min={}, .range_max={}};
  parse_interpolation_methods(&enc, methods);                 
  fprintf(stderr,"Encoding %ix%i image...\n",(int)(img->w),(int)(img->h));
  
  for(int c=0; c<3; c++) {
        fprintf(stderr,"Channel %c: max pixel dist %i; quantization %i",planesYIQ[c],enc.maxD[c],enc.qz[c]/COLOR_STRETCH);
        if (enc.maxD[c] == 0 && enc.qz[c] == 1) fprintf(stderr," (i.e., lossless)");
        fprintf(stderr,"\n");
        if (enc.qz[c]>COLOR_STRETCH && enc.maxD[c]*COLOR_STRETCH < 2*enc.qz[c])
                fprintf(stderr,"  Warning: quantization on channel %c is %i while max distance is only %i!\n",planesYIQ[c],enc.qz[c],enc.maxD[c]);
  }
  FILE *f = fopen(out,"wb");
  fwrite(magic,4,1,f);
  fputc((img->w >> 0) & 0xFF,f);
  fputc((img->w >> 8) & 0xFF,f);
  fputc((img->h >> 0) & 0xFF,f);
  fputc((img->h >> 8) & 0xFF,f);
  fputc((enc.qz[0] >> 0) & 0xFF,f);
  fputc((enc.qz[0] >> 8) & 0xFF,f);
  fputc((enc.qz[1] >> 0) & 0xFF,f);
  fputc((enc.qz[1] >> 8) & 0xFF,f);
  fputc((enc.qz[2] >> 0) & 0xFF,f);
  fputc((enc.qz[2] >> 8) & 0xFF,f);
  fputc((cutoff >> 0) & 0xFF,f);
  fputc((cutoff >> 8) & 0xFF,f);
  fprintf(stderr,"Using RAC range [%g-%g] (cutoff %i)\n",cutoff/4096.0,(4096-cutoff)/4096.0,cutoff);

  symb_init_write(&enc.coder,f,cutoff);
  ctx_init(enc.ctx,cutoff);
//  extra_ctx_init(&enc,cutoff);

/*
  for (int c=0; c<3; c++) {
    for (int m=0; m<4; m++) {
    for (int n=0; n<QZ_BVAR_S; n++) {
      fprintf(stderr,"begin split chances %c others=%i bvar=%i: ",planesYIQ[c],m,n);
      for (int s=0; s<QZ_SIZE_S; s++) {
        fprintf(stderr," %g",1.0*enc.ctx->splitCont[s][c][m][n].chance/4096);
      }
      fprintf(stderr,"\n");
    }
  }
  }
  */
  output_interpolation_methods(&enc);
  for (int i=0; i<skipsplit; i++) symb_put_simple_bit(&enc.coder,1,NULL);
  symb_put_simple_bit(&enc.coder,0,NULL);
  fprintf(stderr,"Areas over %i pixels will always be splitted\n",ms);
  for (int i = 0; i < 16; i++) {
    symb_put_simple_bit(&enc.coder,(imagic >> i) & 1,NULL);
  }
  

  image_set_alpha(img,0);
  rgb2yiq(img);
 
  int min[3] = {255,510,510};
  int max[3] = {0,0,0};
  for(int y = 0; y < img->h; y++) {
     for(int x = 0; x < img->w; x++) {
        pixel_t *p = image_pixel(img,x,y);
        for(int c=0; c<3; c++) {
            if (UNSTRETCH(p->d[c]) < min[c]) min[c]=UNSTRETCH(p->d[c]);
            if (UNSTRETCH(p->d[c]) > max[c]) max[c]=UNSTRETCH(p->d[c]);
        }
     }
  }
  
  for(int c=0; c<3; c++) {
        fprintf(stderr,"Channel %c: range=[%i,%i]\n",planesYIQ[c],min[c],max[c]);
        enc.range_min[c]=min[c];
        enc.range_max[c]=max[c];
//        fputc((enc.range_min[c] >> 0) & 0xFF,f);
//        fputc((enc.range_max[c] >> 0) & 0xFF,f);
        output_min_max(&enc,enc.range_min[c],enc.range_max[c],0,(c==0?255:510),c);

  }
/*  fputc((
        ((enc.range_min[1] >> 8) << 0) |
        ((enc.range_max[1] >> 8) << 1) |
        ((enc.range_min[2] >> 8) << 2) |
        ((enc.range_max[1] >> 8) << 3) )   & 0x0F,f);
  */  
  
  uint32_t p= img->w * img->h;
  uint32_t bvar[3]={};
  output_pixel(&enc,0,0,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  output_pixel(&enc,img->w - 1,0,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  output_pixel(&enc,img->w - 1,img->h - 1,NULL,p,0,NULL,NULL,NULL,NULL,bvar);
  output_pixel(&enc,0,img->h - 1,NULL,p,0,NULL,NULL,NULL,NULL,bvar);

  encode_recurse(&enc,0,img->w - 1,0,img->h - 1,0,0,0,0,0,0,0,0);

  for (int i = 0; i < 16; i++) {
    symb_put_simple_bit(&enc.coder,(imagic >> i) & 1,NULL);
  }
  uint32_t crc=image_crc(&rec);
  for (int i = 0; i < 32; i++) {
    symb_put_simple_bit(&enc.coder,(crc >> i) & 1,NULL);
  }

  symb_flush(&enc.coder);
  fclose(f);

  yiq2rgb(&rec);
  image_save("out.pnm",&rec);
  
  fprintf(stderr,"outputted pixels: %i(%.2f%%) %i(%.2f%%) %i(%.2f%%)\n",enc.outputted_pixels[0],100.0 * enc.outputted_pixels[0] / (img->w * img->h),enc.outputted_pixels[1],100.0 * enc.outputted_pixels[1] / (img->w * img->h),enc.outputted_pixels[2],100.0 * enc.outputted_pixels[2] / (img->w * img->h));
  fprintf(stderr,"Splits: Y:%i  I:%i  Q:%i\n",enc.splits[0],enc.splits[1],enc.splits[2]);
  for (int c=0; c<3; c++) {
    fprintf(stderr,"bits %c: %.1f (%.1f split (%g/symb), %.1f interpol (%g/symb), %.1f data (%g/symb))\n",
        planesYIQ[c],
        (enc.ctx->bitsSplit[c]+enc.ctx->bitsInterpol[c]+enc.ctx->bitsPixeldata[c])*log4k_base,
        enc.ctx->bitsSplit[c]*log4k_base,
        enc.ctx->symbSplit[c]>0 ? enc.ctx->bitsSplit[c]/enc.ctx->symbSplit[c]*log4k_base : 0.0,
        enc.ctx->bitsInterpol[c]*log4k_base,
        enc.ctx->symbInterpol[c]>0 ? enc.ctx->bitsInterpol[c]/enc.ctx->symbInterpol[c]*log4k_base : 0.0,
        enc.ctx->bitsPixeldata[c]*log4k_base,
        enc.ctx->symbPixeldata[c]>0 ? enc.ctx->bitsPixeldata[c]/enc.ctx->symbPixeldata[c]*log4k_base : 0.0
      );
  }
/*  fprintf(stderr,"bits global: %.1f split (%g/symb)\n",
        enc.ctx->bitsSplit[3]*log4k_base,
        enc.ctx->symbSplit[3]>0 ? enc.ctx->bitsSplit[3]/enc.ctx->symbSplit[3]*log4k_base : 0.0
        );*/
//  fprintf(stderr,"Points: left=%i right=%i top=%i bottom=%i\n",enc.pLeft,enc.pRight,enc.pTop,enc.pBottom);
  if (enc.nb_interpolation_methods>1) {
    fprintf(stderr,"Interpolation methods used: \n");
        uint64_t sum=0;
        for(int i=0; i<MAX_INTERPOLATIONS+1; i++) {
         sum += enc.method[i];
        }
        for(int i=0; i<MAX_INTERPOLATIONS+1; i++) {
         if (enc.method[i]>0) {
          fprintf(stderr,"%i. ",i);
          fprintf(stderr,"%s",interpolation_name[i]);
          fprintf(stderr,": %.2f%%  (gain:%e)\n",100.0*enc.method[i]/sum, (double)enc.method_gain[i]);
         }
        }
    fprintf(stderr,"\n");
  }

/*
  fprintf(stderr,"Final border variances histogram: \n");
     for(int c=0; c<3; c++) {
        for(int j=0; j<QZ_BVAR; j++) {
                fprintf(stderr,"bvar[%i]=%i : %li times.\n",c,j,stats_bvar[c][j]);
        }
     }
     */
/*
  fprintf(stderr,"Final chances: \n");
     for(int c=0; c<3; c++) {
                for(int j=0; j<30; j++) {
        for(int i=0; i<32; i++) {
                fprintf(stderr,"chance_c_i_j(%i,%i,%i,%g).\n",c,i,j,1.0*enc.ctx->diff[i][c].chs[j].count/enc.ctx->diff[i][c].chs[j].total);
                }
        }
     }
*/  

/*
  for (int c=0; c<3; c++) {
    for (int m=0; m<4; m++) {
    for (int n=0; n<QZ_BVAR_S; n++) {
      fprintf(stderr,"end split chances %c others=%i bvar=%i: ",planesYIQ[c],m,n);
      for (int s=0; s<QZ_SIZE_S; s++) {
        if (enc.ctx->splitCont[c][s][m][n].count) {
//        fprintf(stderr," %g",1.0*enc.ctx->splitCont[c][s][m][n].count/enc.ctx->splitCont[c][s][m][n].total);
        fprintf(stderr," %g",1.0*enc.ctx->splitCont[s][c][m][n].chance/4096);
        } else {
        fprintf(stderr," /");
        }
      }
      fprintf(stderr,"\n");
    }
  }
  }
  */
  symb_show_stats(stderr,&enc.coder);
        
  image_free(&rec);
}

// main is just main, right?
int main(int argc, char **argv) {
  fprintf(stderr,"Using %lu bytes of context data\n",(unsigned long) sizeof(ctx_t));
  char* methods;
  argc--;
  argv++;
  int mode = 0; // 0=encode 1=decode
  if (argc > 0 && strcmp(argv[0],"-e") == 0) {
    mode = 0;
    argc--;
    argv++;
  }
  if (argc > 0 && strcmp(argv[0],"-d") == 0) {
    mode = 1;
    argc--;
    argv++;
  }
  if (argc < 2 && mode == 0) {
    fprintf(stderr,"Usage: jif [-e] [-i skipsplit] [-c cutoff] [-m METHODS] <input> <output> [quality_y] [quality_i] [quality_q] [factor] [qz_y] [qz_i] [qz_q]\n");
    return 1;
  } else if (argc < 2 && mode == 1) {
    fprintf(stderr,"Usage: jif [-d] <input> <output>\n");
    return 1;
  }
  if (mode == 0) {
    int minsplit=0;
    int cutoff=8;
    if (argc > 0 && strcmp(argv[0],"-i") == 0) {
      if (argc < 3) {fprintf(stderr,"Option -i expects a number\n"); return 1; }
      minsplit=(int)strtol(argv[1],NULL,10);
      argc -= 2;
      argv += 2;
    }
    if (argc > 0 && strcmp(argv[0],"-c") == 0) {
      if (argc < 3) {fprintf(stderr,"Option -c expects a number\n"); return 1; }
      cutoff=(int)strtol(argv[1],NULL,10);
      argc -= 2;
      argv += 2;
      if (cutoff<1 || cutoff>2000) { fprintf(stderr,"Warning: cutoff value should be between 1 and 2000, using default (8)\n"); cutoff=8; }
    }
    if (argc > 0 && strcmp(argv[0],"-m") == 0) {
      if (argc < 3) {fprintf(stderr,"Option -m expects a list of methods\n"); return 1;}
      methods = argv[1];
      argc -= 2;
      argv += 2;
    } else {
      methods = "2";
    }
    // NTSC:  Y: 4 MHz  I: 1.3 MHz  Q: 0.4 MHz
    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : 12;
    double qq = argc > 4 ? strtod(argv[4],NULL) : 20;
    double factor = argc > 5 ? strtod(argv[5],NULL) : 5;
    double qzy = argc > 6 ? strtod(argv[6],NULL) : q/2.0-1.0;
    double qzi = argc > 7 ? strtod(argv[7],NULL) : qi/2.0-1.0;
    double qzq = argc > 8 ? strtod(argv[8],NULL) : qq/2.0-1.0;
    fprintf(stderr,"Reading image file: %s\n",argv[0]);
    image_t image;
    image_load(argv[0],&image);
    encode(&image,q,qi,qq,factor,qzy,qzi,qzq,argv[1],methods,minsplit,cutoff);
    image_free(&image);
//    fprintf(stderr,"maxdists= Y:%.17g  I:%.17g  Q:%.17g\n",maxdist[0],maxdist[1],maxdist[2]);
  } else {
    image_t image;
    decode(&image,argv[0]);
    image_save(argv[1],&image);
    image_free(&image);
  }
  return 0;
}
