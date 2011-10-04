/*  JPIF - an experimental segment-based image format encoder/decoder
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
#include "indexing.h"
#include "contexts.h"

//#define DEBUGMODE

// magic header bytes
uint8_t const magic[] = "JPIF";

// magic marker at the end (inside rac)
uint16_t const imagic = 0x5E2D;




// Hilbert curve shapes: 12 cases
// bit 1: B
// bit 2: rectangular
// bit 4: reverse splitorder
// bit 8: reverse pixelorder
// bit 16: horizontal split
typedef enum {
  SHAPE_0 = 0,
  SHAPE_2 = 16,
#if HILBERT_CURVE == 1
  SHAPE_0A =   2|8|16,
  SHAPE_0B = 1|2|4|16,
  SHAPE_1 = 4|8|16,
  SHAPE_1A = 2|4,
  SHAPE_1B = 1|2|8,
  SHAPE_2A = 2|8,
  SHAPE_2B = 1|2|4,
  SHAPE_3 = 4|8,
  SHAPE_3A = 2|4|16,
  SHAPE_3B = 1|2|8|16
#endif
} shape_t;
const char * shapestr[32] = { 
  [SHAPE_0] = "SHAPE_0",
  [SHAPE_2] = "SHAPE_2",
#if HILBERT_CURVE == 1
  [SHAPE_1] = "SHAPE_1",
  [SHAPE_3] = "SHAPE_3",
  [SHAPE_0A] = "SHAPE_0A",
  [SHAPE_1A] = "SHAPE_1A",
  [SHAPE_2A] = "SHAPE_2A",
  [SHAPE_3A] = "SHAPE_3A",
  [SHAPE_0B] = "SHAPE_0B",
  [SHAPE_1B] = "SHAPE_1B",
  [SHAPE_2B] = "SHAPE_2B",
  [SHAPE_3B] = "SHAPE_3B"
#endif
};
#define SHAPE_SPLITORDER(X) (!!((X)&4))
#define SHAPE_PIXELORDER(X) (!!((X)&8))
#define SHAPE_SPLITDIR(X) (!!((X)&16))


// initial color (should never become visible)
static const pixel_t purple = {.d= {150,255,500},.a=0};

// statistics
double static maxdist[3] = {0,0,0};
long int static stats_bvar[3][QZ_BVAR] = {};

// plane names
static const char planesRGB[4] = "RGBA";
static const char planesYIQ[4] = "YIQA";

long pixels_at_depth[QZ_DEPTH] = {};

uint16_t path[3] = {STRETCH(0),STRETCH(255),STRETCH(255)};
int pixelcounter = 0;




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

// calculate perceptual distance between an area and its interpolation
// (Pasi Fraenti, 1997)
// --> does not seem to give better results at the moment
/*
void static area_distance_new(coder_t *encode, int x1, int x2, int y1, int y2, int *maxd, double *avgd, uint64_t *ndiff) {
  xwpixel_t tot_d1 = {};
  xwpixel_t tot_d2 = {};
  xwpixel_t tot_d3 = {};
  xwpixel_t tot_d4 = {};
  for (int y = y1; y <= y2; y++) {
    for (int x = x1; x <= x2; x++) {
//      fprintf(stderr,"Computing dist for (%i,%i)\n",x,y);
      xwpixel_t sigma_a = image_pixel_contrast(encode->source,x,y);
      xwpixel_t sigma_b = image_pixel_contrast(encode->rec,x,y);
//      fprintf(stderr,"Contrasts computed\n");
      xwpixel_t d1 = {};
      xwwpixel_add(&d1, &sigma_a, 1, 0);
      xwwpixel_add(&d1, &sigma_b, -1, 0);
      xwpixel_square(&d1, 0);

      xwpixel_max(&sigma_a, 1, 0);
      xwpixel_mul(&d1, 100, 0);
      xwwpixel_div(&d1, &sigma_a, 0);
      xwwpixel_add(&tot_d1, &d1, 1, 0);

//      fprintf(stderr,"d1(%i,%i) = [%lli,%lli,%lli]\n",x,y,d1.wd[0],d1.wd[1],d1.wd[2]);

      xwpixel_t gxa = image_pixel_gx(encode->source,x,y);
      xwpixel_t gxb = image_pixel_gx(encode->rec,x,y);
      xwpixel_t gya = image_pixel_gy(encode->source,x,y);
      xwpixel_t gyb = image_pixel_gy(encode->rec,x,y);

      xwwpixel_add(&gxa, &gxb, -1, 0);
      xwpixel_abs(&gxa, 0);
      xwwpixel_add(&gya, &gyb, -1, 0);
      xwpixel_abs(&gya, 0);
      xwpixel_t d2 = {};
      xwwpixel_add(&d2, &gxa, 50, 0);
      xwwpixel_add(&d2, &gya, 50, 0);
      xwwpixel_div(&d2, &sigma_a, 0);
      xwwpixel_add(&tot_d2, &d2, 1, 0);

      xwpixel_t q_a = image_pixel_nb_distinct(encode->source,x,y);
      xwpixel_t q_b = image_pixel_nb_distinct(encode->rec,x,y);
      xwpixel_t d3 = {};
      xwwpixel_add(&d3, &q_a, 1, 0);
      xwwpixel_add(&d3, &q_b, -1, 0);
      xwpixel_square(&d3, 0);
      xwwpixel_add(&tot_d3, &d3, 100, 0);

      pixel_t *real = image_pixel(encode->source,x,y);
      pixel_t *guess = image_pixel(encode->rec,x,y);
      xwpixel_t d4 = {};
      xwpixel_addu(&d4, real,   1, 0);
      xwpixel_addu(&d4, guess, -1, 0);
      xwpixel_abs(&d4, 0);
      xwwpixel_add(&tot_d4, &d4, 100, 0);

    }
  }
  int64_t nb_pixels = ((y2 - y1 + 1) * (x2 - x1 + 1));

//  xwpixel_div(&tot_d1, nb_pixels, 0);
//  xwpixel_div(&tot_d2, nb_pixels, 0);
//  xwpixel_div(&tot_d3, nb_pixels, 0);

  xwpixel_t tot_d = {};
  xwwpixel_add(&tot_d, &tot_d1, 1, 0);
  xwwpixel_add(&tot_d, &tot_d2, 1, 0);
  xwwpixel_add(&tot_d, &tot_d3, 1, 0);
  xwwpixel_add(&tot_d, &tot_d4, 1, 0);
  
//  fprintf(stderr,"area (%i,%i)-(%i,%i): D1 = %lli | D2 = %lli | D3 = %lli | D4 = %lli | D = %lli\n",x1,y1,x2,y2,tot_d1.wd[0],tot_d2.wd[0],tot_d3.wd[0],tot_d4.wd[0],tot_d.wd[0]);

  maxd[0] = tot_d.wd[0];
  maxd[1] = tot_d.wd[0];
  maxd[2] = tot_d.wd[0];

  ndiff[0] = tot_d.wd[0];
  ndiff[1] = tot_d.wd[1];
  ndiff[2] = tot_d.wd[2];
}
*/

// calculate simple pixel-by-pixel color distance between an area and its interpolation
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
      fd[0] += distance[0];
      fd[1] += distance[1];
      fd[2] += distance[2];
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


//                  A-----X-----B
//                  |           |
//                  Y           |
//                  |           |
//                  C-----?-----D
void guess_median(pixel_t *guess, pixel_t *A, pixel_t *B, pixel_t *C, pixel_t *D, pixel_t *X, pixel_t *Y, pixel_t *CD, int mask, coder_t *coder) {
        wpixel_t XC_A = {};
        wpixel_add(&XC_A, X, 1, mask);
        wpixel_add(&XC_A, C, 1, mask);
        wpixel_add(&XC_A, A,-1, mask);
        wpixel_t XD_B = {};
        wpixel_add(&XD_B, X, 1, mask);
        wpixel_add(&XD_B, D, 1, mask);
        wpixel_add(&XD_B, B,-1, mask);
        pixel_t XCA = {};
        pixel_t XDB = {};
        wpixel_set(&XCA,&XC_A,mask);
        wpixel_set(&XDB,&XD_B,mask);
//        round_color(&coder->color_data,&XCA);
//        round_color(&coder->color_data,&XDB);
//        round_color(&coder->color_data,CD);
//        pixel_median7(guess,CD,&XCA,&XDB,C,D,A,B,mask);
//        pixel_median7(guess,CD,&XCA,&XDB,C,D,X,Y,mask);
        pixel_median5(guess,CD,&XCA,&XDB,C,D,mask);
//        pixel_median3(guess,CD,&XCA,&XDB,mask);

        round_color(&coder->color_data,guess);
}

// order of guess: left, right, top, bottom
void inline static border_guess(int guess_method, coder_t *coder, int px1, int px2, int py1, int py2, int xm, int ym, pixel_t *guess, int mask, int corner) {

          if (corner == 0) {    
            pixel_linear(&guess[0],image_pixel(coder->rec,px1,py1),ym - py1,image_pixel(coder->rec,px1,py2),py2 - ym,mask);
            return;
          }
          if (corner == 2) {
            pixel_linear(&guess[2],image_pixel(coder->rec,px1,py1),xm - px1,image_pixel(coder->rec,px2,py1),px2 - xm,mask);
            return;
          }

          pixel_t * A = image_pixel(coder->rec,px1,py1);
          pixel_t * B = image_pixel(coder->rec,px2,py1);
          pixel_t * C = image_pixel(coder->rec,px1,py2);
          pixel_t * D = image_pixel(coder->rec,px2,py2);

          if (corner == 3) {
            pixel_linear(&guess[2],image_pixel(coder->rec,px1,py1),xm - px1,image_pixel(coder->rec,px2,py1),px2 - xm,mask);
            pixel_t * X;
            if (py1 == 0) {X = &guess[2];} else {X = image_pixel(coder->rec,xm,py1);}
            pixel_t CD = {};
            pixel_linear(&CD,image_pixel(coder->rec,px1,py2),xm - px1,image_pixel(coder->rec,px2,py2),px2 - xm,mask);
  //          guess_median(&guess[3],A,B,C,D,X,Y,&CD,mask,coder);
            guess_median(&guess[3],A,B,C,D,X,NULL,&CD,mask,coder);
            return;
          }
          if (corner == 1) {
            pixel_linear(&guess[0],image_pixel(coder->rec,px1,py1),ym - py1,image_pixel(coder->rec,px1,py2),py2 - ym,mask);
            pixel_t * Y;
            if (px1 == 0) {Y = &guess[0];} else {Y = image_pixel(coder->rec,px1,ym);}
            pixel_t BD = {};
            pixel_linear(&BD,image_pixel(coder->rec,px2,py1),ym - py1,image_pixel(coder->rec,px2,py2),py2 - ym,mask);
//          guess_median(&guess[1],A,C,B,D,Y,X,&BD,mask,coder);
            guess_median(&guess[1],A,C,B,D,Y,NULL,&BD,mask,coder);
            return;
          }

/* old stuff:
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
          pixel_set(&guess[0],image_pixel(coder->rec,px1,ym),mask);
          pixel_set(&guess[1],image_pixel(coder->rec,px2,ym),mask);
          pixel_set(&guess[2],image_pixel(coder->rec,xm,py1),mask);
          pixel_set(&guess[3],image_pixel(coder->rec,xm,py2),mask);
          return;
  }
*/  
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

// don't take sides into account of length X or longer (too expensive)
void static compute_bvar(coder_t *encode, int x1, int y1, int x2, int y2, int mask, uint32_t *bvar, int hastop, int hasright, int hasbottom, int hasleft) { //, pixel_t *avg) {
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  int val;

  for(int c=0;c<3;c++) {
   if(!(mask & 1<<c)) {
    int64_t sum=0;
    int64_t sumq=0;
    int64_t count=0;
    if (hasleft) {
     int y_step = (py2-y1)/50 + 1;
     for(int y=y1+y_step/2; y<py2; y += y_step) {
         val=UNSTRETCH(image_pixel(encode->rec,px1,y)->d[c]);
         sum += val;
         sumq += val*val;
         count++;
     }
    }
    if (hastop) {
     int x_step = (px2-x1)/50 + 1;
     for(int x=x1+x_step/2; x<px2; x += x_step) {
         val=UNSTRETCH(image_pixel(encode->rec,x,py1)->d[c]);
         sum += val;
         sumq += val*val;
         count++;
     }
    }

    val=UNSTRETCH(image_pixel(encode->rec,px1,py1)->d[c]);
    sum += val;
    sumq += val*val;
    val=UNSTRETCH(image_pixel(encode->rec,px1,py2)->d[c]);
    sum += val;
    sumq += val*val;
    val=UNSTRETCH(image_pixel(encode->rec,px2,py1)->d[c]);
    sum += val;
    sumq += val*val;
    val=UNSTRETCH(image_pixel(encode->rec,px2,py2)->d[c]);
    sum += val;
    sumq += val*val;
    count += 4;
    
    if (count>0) {
      int scale = 1<<14;
      if (c==0) scale = 1<<16;
      uint32_t var=((sumq*scale)-(sum*sum*scale/count))/count;
      bvar[c] = var;
      assert(bvar[c]>=0);
//      if (bvar[c] < 0) bvar[c]=0;
    } else {
      bvar[c] = 1<<16;
    }
//    stats_bvar[c][qbvar(bvar[c],QZ_BVAR)]++;
   }
  }
}

void static compute_bdiffE(coder_t *encode, int x1, int y1, int x2, int y2, int mask, uint32_t *bdiffE, int hastop, int hasright, int hasbottom, int hasleft) {
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  pixel_t interpol = {};
  bdiffE[0] = 0;
  bdiffE[1] = 0;
  bdiffE[2] = 0;

  if (hasleft) {
     int y_step = (py2-y1)/50 + 1;
     for(int y=y1+y_step/2; y<py2; y += y_step) {
         pixel_linear(&interpol,image_pixel(encode->rec,px1,py1),y - py1,image_pixel(encode->rec,px1,py2),py2 - y,mask); 
         for (int c=0; c<3; c++) bdiffE[c] += UNSTRETCH(abs(image_pixel(encode->rec,px1,y)->d[c] - interpol.d[c]));
     }
  } else {
      for (int c=0; c<3; c++) bdiffE[c] = 50*(c==0?50:100);
  }
  if (hastop) {
     int x_step = (px2-x1)/50 + 1;
     for(int x=x1+x_step/2; x<px2; x += x_step) {
         pixel_linear(&interpol,image_pixel(encode->rec,px1,py1),x - px1,image_pixel(encode->rec,px2,py1),px2 - x,mask); 
         for (int c=0; c<3; c++) bdiffE[c] += UNSTRETCH(abs(image_pixel(encode->rec,x,py1)->d[c] - interpol.d[c]));
     }
  } else {
     for (int c=0; c<3; c++) bdiffE[c] = 50*(c==0?50:100);
  }
}



#if FULL_BORDERS
// recursive function for encoding
void static encode_recurse_line(coder_t *encode, int x1, int x2, int y1, int y2, int mask, int depth, uint32_t *bvar, int shape) {
  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  uint32_t p = w * h;
  int size = ilog2(p);
  int qsize_sl = quantize_log_uint32(p,QZ_SIZE_SL+1)-1;
  int current_method = 0;
  assert(x2 >= x1 && y2 >= y1);
  assert(x2 == x1 || y2 == y1);
  
  if (p < 3) return;

//  fprintf(stderr,"Line from (%i,%i) to (%i,%i) at depth %i\n",x1,y1,x2,y2,depth);
  pixel_t *orig[4] = {
        image_pixel(encode->rec,x1,y1),               //left top
        image_pixel(encode->rec,x2,y1),               //right top
        image_pixel(encode->rec,x1,y2),               //left bottom
        image_pixel(encode->rec,x2,y2)                //right bottom
        };
  int always_split[3] = {0,0,0};
  for (int c=0; c<3; c++) {
     if (!(mask & (1 << c))) {
       if (abs(orig[0]->d[c]-orig[3]->d[c]) > ALWAYS_SPLIT_LINE_DIST*encode->qz[c]) {
                mask |= 1<<c;
                always_split[c] = 1;
        }                
     }
  }
  
  int method;
  if (p-2 > encode->minsplitL) {
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
      interpolation(encode->interpolation_methods[i],encode->rec,x1,x2,y1,y2,mask,0,0,0,0,&encode->color_data);
      if (x1==x2) area_distance(encode,x1,x1,y1+1,y2-1,maxd,&avgd,rd);
      if (y1==y2) area_distance(encode,x1+1,x2-1,y1,y1,maxd,&avgd,rd);
//      area_distance(encode,x1,x2,y1,y2,maxd,&avgd,rd);
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

    if (factor > MAX_AVG_ERROR_LINE*p) factor = MAX_AVG_ERROR_LINE*p;

    double f = 1.0;
    if (p>3) f = LINE_FACTOR;
    // stop criterion
    for (int c = 0; c < 3 ; c++) {
      if (!(mask & (1 << c))) {
        method=maxd_i[c];
//        fprintf(stderr,"[%i] newmask %i   others %i\n",c,newmask,others);
        if (   (f*maxd[c] <= encode->maxD[c] )   || rd[c] <= encode->maxD[c] * factor ) {
          if (maxd[c] > encode->maxD[c]) method=rd_i[c];
          if (maxd[c] > maxdist[c]) maxdist[c] = maxd[c];
          interpolation(encode->interpolation_methods[method],encode->rec,x1,x2,y1,y2,~(1<<c),0,0,0,0,&encode->color_data);
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

    output_mask(encode,mask,newmask,qsize_sl,depth,bvar,1,orig[0],orig[3]);

     for (int c = 0; c < 3 ; c++) {
      if (!(mask & (1 << c))) {
        method=maxd_i[c];
        if ((mask & (1 << c)) != (newmask & (1 << c))) {
          update_avg_split_depths(encode,c,depth,1);
          output_interpolation_choice(encode,size,c,method,current_method);
        }
      }
     }
    

    mask = newmask;
  }

  for (int c=0; c<3; c++) {
     if (always_split[c]) mask &= ~(1<<c);
  }

  if (mask == 7) return;

  pixel_t guess[4];

  int xm = (x1 + x2) / 2 ;
  int ym = (y1 + y2) / 2 ;
  assert(xm >= 0);
  assert(ym >= 0);


  depth++;
  assert(depth < DEPTH_LEVELS);

  
  if(x1==x2) {
  // horizontal
        border_guess(current_method,encode,x1,x2,y1,y2,xm,ym,guess,0,1);
        output_pixel(encode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar,depth,0,0);
        if(SHAPE_PIXELORDER(shape)) {
          //bottom-top
          encode_recurse_line(encode,x1,x1,ym,y2,mask,depth,bvar,shape); 
          encode_recurse_line(encode,x1,x1,y1,ym,mask,depth,bvar,shape); 
        } else {
          //top-bottom
          encode_recurse_line(encode,x1,x1,y1,ym,mask,depth,bvar,shape); 
          encode_recurse_line(encode,x1,x1,ym,y2,mask,depth,bvar,shape); 
        }
  } else {
  // vertical
        assert(y1==y2);
        border_guess(current_method,encode,x1,x2,y1,y2,xm,ym,guess,0,3);
        output_pixel(encode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar,depth,0,1);
        if(SHAPE_PIXELORDER(shape)) {
          //right-left
          encode_recurse_line(encode,xm,x2,y1,y1,mask,depth,bvar,shape); 
          encode_recurse_line(encode,x1,xm,y1,y1,mask,depth,bvar,shape); 
        } else {
          //left-right
          encode_recurse_line(encode,x1,xm,y1,y1,mask,depth,bvar,shape); 
          encode_recurse_line(encode,xm,x2,y1,y1,mask,depth,bvar,shape); 
        }
  }
}
#endif

#if FULL_BORDERS == 0
int static inline nn(int x,int y,int px1,int px2,int py1,int py2) {
        if (x == px1 && y == py1) return 0;
        if (x == px1 && y == py2) return 0;
        if (x == px2 && y == py1) return 0;
        if (x == px2 && y == py2) return 0;
        return 1;
}
int static inline nn2(int x,int y,int px1,int px2,int x1,int py1,int py2,int y1,int other) {
        if (x == px1 && y == py1) return 0;
        if (x == px1 && y == py2) return 0;
        if (x == px2 && y == py1) return 0;
        if (x == px2 && y == py2) return 0;
        if (x == x1  && y == y1 && other) return 0;
        return 1;
}
#endif

#if FFV1
void static encode_ffv1(coder_t *encode, int x1, int x2, int y1, int y2) {
  pixel_t black = {.d= {0,0,0}};
  for (int y = y1; y <= y2; y++) {
    if (y % (y2/((1<<PROGRESS_BAR_VERBOSITY)-1)) == 0) fprintf(stderr,PROGRESS_BAR_CHAR);

    for (int x = x1; x <= x2; x++) {
        pixel_t *orig[6] = {};
        if (x>0 && y>0) {
           orig[0] = image_pixel(encode->rec,x-1,y-1);  //top left
        } else {
           orig[0] = &black;
        }
        if (x>0) {
           orig[1] = image_pixel(encode->rec,x-1,y);    //left
        } else {
           if (y>0) {
             orig[1] = image_pixel(encode->rec,x,y-1); 
           } else {
             orig[1] = &black;
           }
        }
        if (y>0) {
           orig[2] = image_pixel(encode->rec,x,y-1);    //top
        } else {
           orig[2] = &black;
        }
        if (x>1) {
           orig[3] = image_pixel(encode->rec,x-2,y);    //left left
        } else {
           if (x>0 && y>0) {
             orig[3] = image_pixel(encode->rec,x-1,y-1); 
           } else {
             orig[3] = &black;
           }
        }
        if (y>1) {
           orig[4] = image_pixel(encode->rec,x,y-2);    //top top
        } else {
           orig[4] = &black;
        }
        if (y>0 && x<x2) {
           orig[5] = image_pixel(encode->rec,x+1,y-1);    //top right
        } else {
           orig[5] = orig[2];
        }
        wpixel_t myguess = {};
        wpixel_add(&myguess,orig[1],1,0);
        wpixel_add(&myguess,orig[2],1,0);
        wpixel_add(&myguess,orig[0],-1,0);
        pixel_t grad;
        wpixel_set(&grad,&myguess,0);
        pixel_t *pixels[3] = {orig[1],orig[2],&grad};
        pixel_t guess;
        for (int c=0; c<3; c++) {
           guess.d[c] = median_3(pixels[0]->d[c],pixels[1]->d[c],pixels[2]->d[c]);
        }
//        round_color(&encode->color_data,&guess);
        output_pixel_ffv1(encode,x,y,&guess,0,orig[0],orig[1],orig[2],orig[3],orig[4],orig[5]);
    }     
  }
}
void static decode_ffv1(coder_t *decode, int x1, int x2, int y1, int y2) {
  pixel_t black = {.d= {0,0,0}};
  for (int y = y1; y <= y2; y++) {
    if (y % (y2/((1<<PROGRESS_BAR_VERBOSITY)-1)) == 0) fprintf(stderr,PROGRESS_BAR_CHAR);
    for (int x = x1; x <= x2; x++) {
        pixel_t *orig[6] = {};
        if (x>0 && y>0) {
           orig[0] = image_pixel(decode->rec,x-1,y-1);  //top left
        } else {
           orig[0] = &black;
        }
        if (x>0) {
           orig[1] = image_pixel(decode->rec,x-1,y);    //left
        } else {
           if (y>0) {
             orig[1] = image_pixel(decode->rec,x,y-1); 
           } else {
             orig[1] = &black;
           }
        }
        if (y>0) {
           orig[2] = image_pixel(decode->rec,x,y-1);    //top
        } else {
           orig[2] = &black;
        }
        if (x>1) {
           orig[3] = image_pixel(decode->rec,x-2,y);    //left left
        } else {
           if (x>0 && y>0) {
             orig[3] = image_pixel(decode->rec,x-1,y-1); 
           } else {
             orig[3] = &black;
           }
        }
        if (y>1) {
           orig[4] = image_pixel(decode->rec,x,y-2);    //top top
        } else {
           orig[4] = &black;
        }
        if (y>0 && x<x2) {
           orig[5] = image_pixel(decode->rec,x+1,y-1);    //top right
        } else {
           orig[5] = orig[2];
        }
        wpixel_t myguess = {};
        wpixel_add(&myguess,orig[1],1,0);
        wpixel_add(&myguess,orig[2],1,0);
        wpixel_add(&myguess,orig[0],-1,0);
        pixel_t grad;
        wpixel_set(&grad,&myguess,0);
        pixel_t *pixels[3] = {orig[1],orig[2],&grad};
        pixel_t guess;
        for (int c=0; c<3; c++) {
           guess.d[c] = median_3(pixels[0]->d[c],pixels[1]->d[c],pixels[2]->d[c]);
        }
        input_pixel_ffv1(decode,x,y,&guess,0,orig[0],orig[1],orig[2],orig[3],orig[4],orig[5]);
    }     
  }
}

#endif

static inline int pixel_compute_dist(pixel_t *a, pixel_t *b, int mask) {
       wpixel_t dist = {};
       wpixel_addu(&dist,a,1,mask);
       wpixel_addu(&dist,b,-1,mask);
       wpixel_abs(&dist,mask);
       int total=0;
       for (int c=0; c<3; c++) total+= dist.wd[c];
       return total;
}

static inline int pixel_dist(image_t *img, int xa, int ya, int xb, int yb, int mask) {
       return pixel_compute_dist(image_pixel(img,xa,ya),image_pixel(img,xb,yb),mask);
}

#define MAX_WIDTH_HEIGHT_RATIO 8
#define MIN_SIZE_DO_HEURISTIC 100

#define KNOWN_MIDDLE_PIXEL_DIFF_WEIGHT 5
#define MAX_KNOWN_BORDER_DIFF_WEIGHT 2
#define MIN_DIFF_DIFF 30
#define ASYM_SPLIT_THRESHOLD 25
#define MIN_BORDER_SIZE_ASYM_SPLIT_CHECK 8
#define MAX_BORDER_SIZE_ASYM_SPLIT_CHECK 50

// asymmetric splitting does not seem to be a very good idea
#define ASYM_SPLITTING 0
//#define ASYM_SPLIT_DO_THRESHOLD 320
//#define ASYM_SPLIT_DO_THRESHOLD 350
#define ASYM_SPLIT_DO_THRESHOLD 250
int asym_splits = 0;

inline static int decide_splitdir(image_t *img,int x1,int x2,int y1,int y2,int px1, int px2, int py1, int py2, int mask,int *xm,int *ym) {
/*  if (x2-x1 > y2-y1) {
        return SHAPE_0;
  } else {
        return SHAPE_2;
  }
  */

  if (x2==x1) return SHAPE_2;
  if (y2==y1) return SHAPE_0;
  if (x2-x1 > MAX_WIDTH_HEIGHT_RATIO*(y2-y1)) return SHAPE_0;
  if (y2-y1 > MAX_WIDTH_HEIGHT_RATIO*(x2-x1)) return SHAPE_2;
  if ((x2-x1+1)*(y2-y1+1) > MIN_SIZE_DO_HEURISTIC) {

  int lrt = pixel_dist(img, px1, py1, px2, py1, mask);  // left-right diff (top)
  int lrb = pixel_dist(img, px1, py2, px2, py2, mask);  // left-right diff (bottom)
  int tbl = pixel_dist(img, px1, py1, px1, py2, mask);  // top-bottom diff (left)
  int tbr = pixel_dist(img, px2, py1, px2, py2, mask);  // top-bottom diff (right)
  
  int lre = 0; // left-right diff with expected (interpolation);
  int lr_max = ASYM_SPLIT_THRESHOLD;
  int lr_min = -1;
  int lr_best = *xm;
  if (py1<y1) {  // top is known
        pixel_t interpol = {};
        pixel_linear(&interpol,image_pixel(img,px1,py1),*xm - px1,image_pixel(img,px2,py1),px2 - *xm,mask);
        lre = KNOWN_MIDDLE_PIXEL_DIFF_WEIGHT*pixel_compute_dist(&interpol, image_pixel(img,*xm,py1), mask);
        if (x2-x1 > MIN_BORDER_SIZE_ASYM_SPLIT_CHECK && x2-x1 < MAX_BORDER_SIZE_ASYM_SPLIT_CHECK) {
         for (int x=px1; x < px2; x++) {
           int diff = pixel_dist(img, x,py1, x+1, py1, mask);
           if (lr_min == -1 || diff < lr_min) lr_min=diff;
           if (diff > lr_max) {
                lr_max = diff;
                lr_best = x+1;
//                if (x<*xm) lr_best++;
           }
         } 
        }
  }
  int tbe = 0; // top-bottom diff with expected (interpolation);
  int tb_max = ASYM_SPLIT_THRESHOLD;
  int tb_min = -1;
  int tb_best = *ym;
  if (px1<x1) {  // left is known
        pixel_t interpol = {};
        pixel_linear(&interpol,image_pixel(img,px1,py1),*ym - py1,image_pixel(img,px1,py2),py2 - *ym,mask);
        tbe = KNOWN_MIDDLE_PIXEL_DIFF_WEIGHT*pixel_compute_dist(&interpol, image_pixel(img,px1,*ym), mask);
        if (y2-y1 > MIN_BORDER_SIZE_ASYM_SPLIT_CHECK && y2-y1 < MAX_BORDER_SIZE_ASYM_SPLIT_CHECK) {
         for (int y=py1; y < py2; y++) {
           int diff = pixel_dist(img, px1,y, px1, y+1, mask);
           if (tb_min == -1 || diff < tb_min) tb_min=diff;
           if (diff > tb_max) {
                tb_max = diff;
                tb_best = y+1;
//                if (y<*ym) tb_best++;
           }
         }
        }
  }

  lr_max -= lr_min;
  tb_max -= tb_min;

#if ASYM_SPLITTING
//  if (lr_max > ASYM_SPLIT_DO_THRESHOLD) { *xm = clamp(lr_best,x1+(x2-x1)/16,x2-(x2-x1)/16); asym_splits++; }
//  if (tb_max > ASYM_SPLIT_DO_THRESHOLD) { *ym = clamp(tb_best,y1+(y2-y1)/16,y2-(y2-y1)/16); asym_splits++; }
  if (lr_max > ASYM_SPLIT_DO_THRESHOLD) { *xm = clamp(lr_best,x1+1+(x2-x1)/16,x2-1-(x2-x1)/16); asym_splits++; }
  if (tb_max > ASYM_SPLIT_DO_THRESHOLD) { *ym = clamp(tb_best,y1+1+(x2-x1)/16,y2-1-(x2-x1)/16); asym_splits++; }
#endif

  lr_max *= MAX_KNOWN_BORDER_DIFF_WEIGHT;
  tb_max *= MAX_KNOWN_BORDER_DIFF_WEIGHT;
  
  if (lrt+lrb+lre+lr_max > MIN_DIFF_DIFF+tbl+tbr+tbe+tb_max) return SHAPE_0;
  if (tbl+tbr+tbe+tb_max > MIN_DIFF_DIFF+lrt+lrb+lre+lr_max) return SHAPE_2;

  if (px1<x1 && py1==y1) return SHAPE_0;        // left known but not top
  if (px1==x1 && py1<y1) return SHAPE_2;        // top known but not bottom
  }
  
  // no strong reason to split in one direction or the other, so use default (keep rectangles square)
  if (x2-x1 > y2-y1) {
        return SHAPE_0;
  } else {
        return SHAPE_2;
  }

/*  
  if (x2-x1 > y2-y1) {          //                      |
        shape = SHAPE_0;        // vertical split     L | R
                                //                      |
                                
  } else {                      //                      T
        shape = SHAPE_2;        // horizontal split   -----
  }                             //                      B
*/

}


// recursive function for encoding
void static encode_recurse(coder_t *encode, int x1, int x2, int y1, int y2, int mask, int previous_method, int hastop, int hasright, int hasbottom, int hasleft, shape_t shape, int depth, int checksplit) {


#if FULL_BORDERS
  int px1 = x1 - 1;
  int py1 = y1 - 1;
  int px2 = x2 + 1;
  int py2 = y2 + 1;
  int myhastop=1;
  int myhasright=1;
  int myhasbottom=1;
  int myhasleft=1;
#else
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  int myhastop=hastop;
  int myhasright=hasright;
  int myhasbottom=hasbottom;
  int myhasleft=hasleft;
#endif  
  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  uint32_t p = w * h;


  if (depth == PROGRESS_BAR_VERBOSITY) fprintf(stderr,PROGRESS_BAR_CHAR);



#ifdef DEBUGMODE
  fprintf(stderr,"%.*sarea (%i,%i)-(%i,%i) ",depth*2,"                                                                               ",x1,y1,x2,y2);
  fprintf(stderr,"shape:%s ",shapestr[shape]);
  if (hastop) fprintf(stderr,"top ");
  if (hasright) fprintf(stderr,"right ");
  if (hasbottom) fprintf(stderr,"bottom ");
  if (hasleft) fprintf(stderr,"left ");
  fprintf(stderr,"\n");
#endif


  if (px2 - px1 <= 1 && py2 - py1 <= 1) return;
  if (x1 == x2+1 || y1 == y2+1) return;

//  fprintf(stderr,"XXXX x=%i y=%i  px=%i py=%i  depth=%i  shape=%s \n",x2-x1+1,y2-y1+1,px2-px1+1,py2-py1+1,depth,shapestr[shape]);

  assert(x2 >= x1 && y2 >= y1);

  uint32_t bvar[3] = {};
//  pixel_t firstguess;
  compute_bvar(encode,x1,y1,x2,y2,0,&bvar[0],myhastop,myhasright,myhasbottom,myhasleft); //,&firstguess);

  uint32_t bdiffE[3] = {};
  compute_bdiffE(encode,x1,y1,x2,y2,0,&bdiffE[0],myhastop,myhasright,myhasbottom,myhasleft); 




  pixel_t *orig[4] = {
        image_pixel(encode->rec,px1,py1),               //left top
        image_pixel(encode->rec,px2,py1),               //right top
        image_pixel(encode->rec,px1,py2),               //left bottom
        image_pixel(encode->rec,px2,py2)                //right bottom
        };





  int method;
  int current_method = previous_method;

  

#if ALWAYS_SPLIT_HEURISTIC
  int always_split[3] = {0,0,0};
  for (int c=0; c<3; c++) {
     if (!(mask & (1 << c))) {
       if (bvar[c] > ALWAYS_SPLIT_AREA_DIST*encode->qz[c]) {
                mask |= 1<<c;
                always_split[c] = 1;
        }                
     }
  }
#endif


  if (p > encode->minsplitA && checksplit && depth%SPLIT_MODULO == 0) {
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
      interpolation(encode->interpolation_methods[i],encode->rec,x1,x2,y1,y2,mask,myhastop,myhasright,myhasbottom,myhasleft,&encode->color_data);
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
    if (factor > p*encode->factor_s) factor = p*encode->factor_s;

    // uncomment to have average error factor
    // factor *= p;

    // stop criterion
    for (int c = 0; c < 3 ; c++) {
      if (!(mask & (1 << c))) {
        method=maxd_i[c];
//        fprintf(stderr,"[%i] newmask %i   others %i\n",c,newmask,others);
        if (   (maxd[c] <= encode->maxD[c] )   || rd[c] <= encode->maxD[c] * factor ) {
//        if (   rd[c] <= encode->maxD[c] * factor ) {
          if (maxd[c] > encode->maxD[c]) method=rd_i[c];
          if (maxd[c] > maxdist[c]) maxdist[c] = maxd[c];
//          interpolation(encode->interpolation_methods[method],encode->rec,x1,x2,y1,y2,~(1<<c),myhastop,myhasright,myhasbottom,myhasleft,&encode->color_data);
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

#if SPLIT_CHANNELS_SEPARATELY == 0
    if (newmask < 7) newmask = 0;
#endif


  int qsize = quantize_log_uint32(p,QZ_SIZE);
//  int qsize_s = quantize_log_uint32(p,QZ_SIZE_S);
    output_mask(encode,mask,newmask,qsize,depth,bvar,bdiffE,0,NULL,NULL);

     for (int c = 0; c < 3 ; c++) {
      if (!(mask & (1 << c))) {
        method=maxd_i[c];
        if ((mask & (1 << c)) != (newmask & (1 << c))) {
          update_avg_split_depths(encode,c,depth,0);
          output_interpolation_choice(encode,qsize,c,method,current_method);
        }
      }
     }
    

    mask = newmask;
  }

#if ALWAYS_SPLIT_HEURISTIC
  for (int c=0; c<3; c++) {
     if (always_split[c]) mask &= ~(1<<c);
  }
#endif


//  if (px1==29 && py1==29) return;

  if (mask == 7) {
//   for (int x = px1; x<x2; x++) pixel_set(image_pixel(encode->source,x,py1),&purple,0);
//   for (int y = py1; y<y2; y++) pixel_set(image_pixel(encode->source,px1,y),&purple,0);
   return;
  }

  int xm = (x1 + x2) / 2 ;
  int ym = (y1 + y2) / 2 ;

#if HILBERT_CURVE == 0
  shape = decide_splitdir(encode->rec,x1,x2,y1,y2,px1,px2,py1,py2,mask,&xm,&ym); 
/*
  if (x2-x1 > y2-y1) {          //                      |
        shape = SHAPE_0;        // vertical split     L | R
                                //                      |

  } else {                      //                      T
        shape = SHAPE_2;        // horizontal split   -----
  }                             //                      B
*/
#endif 


#if FULL_BORDERS
  if ( x1==x2 && y1==y2) {
    // one pixel

    wpixel_t myguess = {};
    wpixel_add(&myguess,image_pixel(encode->rec,x1,py1),2,0);
    wpixel_add(&myguess,image_pixel(encode->rec,x1,py2),2,0);
    wpixel_add(&myguess,image_pixel(encode->rec,px1,y1),2,0);
    wpixel_add(&myguess,image_pixel(encode->rec,px2,y1),2,0);
    wpixel_add(&myguess,image_pixel(encode->rec,px1,py1),-1,0);
    wpixel_add(&myguess,image_pixel(encode->rec,px2,py1),-1,0);
    wpixel_add(&myguess,image_pixel(encode->rec,px1,py2),-1,0);
    wpixel_add(&myguess,image_pixel(encode->rec,px2,py2),-1,0);
    wpixel_div(&myguess,4,0);
    pixel_t grad;
    wpixel_set(&grad,&myguess,0);
    pixel_t *pixels[5] = {
            image_pixel(encode->rec,x1,py1),
            image_pixel(encode->rec,x1,py2),
            image_pixel(encode->rec,px1,y1),
            image_pixel(encode->rec,px2,y1),
            &grad };
    pixel_t guess;
//    int weights[5] = {1,1,1,1,1};
//    pixel_med(&guess, pixels, weights, 5, mask);
    for (int c=0; c<3; c++) {
        guess.d[c] = median5(pixels[0]->d[c],pixels[1]->d[c],pixels[2]->d[c],pixels[3]->d[c],pixels[4]->d[c]);
    }
    output_pixel(encode,x1,y1,&guess,p,mask,orig[0],orig[1],orig[2],orig[3],bvar,depth,1,0);
    return;
  }
  
  int use_firstguess = 1;
  if (w < h-2 || w > h+2) use_firstguess = 0;  
  // if (w == 1 && h == 1) use_firstguess = 1;

  if (use_firstguess && USE_BORDER_GUESSING) {
  
//  fprintf(stderr,"Getting middle pixel %i,%i\n",xm,ym);
  output_pixel(encode,xm,ym,&firstguess,p,mask,orig[2],orig[1],orig[0],orig[3],bvar,bdiffE,depth,1,0);
  if(SHAPE_SPLITDIR(shape)) {
  // horizontal
        if(SHAPE_PIXELORDER(shape)) {
          //right-left
          encode_recurse_line(encode,xm,px2,ym,ym,mask,depth,bvar,shape); 
          encode_recurse_line(encode,px1,xm,ym,ym,mask,depth,bvar,shape); 
        } else {
          //left-right
          encode_recurse_line(encode,px1,xm,ym,ym,mask,depth,bvar,shape); 
          encode_recurse_line(encode,xm,px2,ym,ym,mask,depth,bvar,shape); 
        }
  } else {
  // vertical
        if(SHAPE_PIXELORDER(shape)) {
          //bottom-top
          encode_recurse_line(encode,xm,xm,ym,py2,mask,depth,bvar,shape); 
          encode_recurse_line(encode,xm,xm,py1,ym,mask,depth,bvar,shape); 
        } else {
          //top-bottom
          encode_recurse_line(encode,xm,xm,py1,ym,mask,depth,bvar,shape); 
          encode_recurse_line(encode,xm,xm,ym,py2,mask,depth,bvar,shape); 
        }
  }
  
  } else {
  if(SHAPE_SPLITDIR(shape)) {
  // horizontal
          encode_recurse_line(encode,px1,px2,ym,ym,mask,depth,bvar,shape); 
  } else {
  // vertical
          encode_recurse_line(encode,xm,xm,py1,py2,mask,depth,bvar,shape); 
  }
  }
  depth++;
  checksplit = 1;
  
  switch(shape) {
        case SHAPE_0:
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft, SHAPE_0A,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0B,depth,checksplit); // right side
          break;
        case SHAPE_1:
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,0,hasright,hasbottom,hasleft,SHAPE_1A,depth,checksplit); // bottom side
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   SHAPE_1B,depth,checksplit); // top side
          break;
        case SHAPE_2:
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,0,hasleft,   SHAPE_2A,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2B,depth,checksplit); // bottom side
          break;
        case SHAPE_3:
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, SHAPE_3A,depth,checksplit); // right side
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,  SHAPE_3B,depth,checksplit); // left side
          break;
        case SHAPE_0A:
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,0,hasleft,    SHAPE_2,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft, SHAPE_0,depth,checksplit); // bottom side
          break;
        case SHAPE_0B:
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,0,hasright,hasbottom,hasleft, SHAPE_0,depth,checksplit); // bottom side
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,    SHAPE_1,depth,checksplit); // top side
          break;
        case SHAPE_1A:
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,  SHAPE_3,depth,checksplit); // right side
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          break;
        case SHAPE_1B:
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,  SHAPE_0,depth,checksplit); // right side
          break;
        case SHAPE_2A:
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_0,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,  SHAPE_2,depth,checksplit); // right side
          break;
        case SHAPE_2B:
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,  SHAPE_2,depth,checksplit); // right side
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,   SHAPE_3,depth,checksplit); // left side
          break;
        case SHAPE_3A:
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,0,hasright,hasbottom,hasleft, SHAPE_1,depth,checksplit); // bottom side
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,    SHAPE_3,depth,checksplit); // top side
          break;
        case SHAPE_3B:
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,0,hasleft,    SHAPE_3,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft, SHAPE_2,depth,checksplit); // bottom side
          break;
  }
#else

#if HILBERT_CURVE == 1
  xm = (px1+px2)/2;
  ym = (py1+py2)/2;
#endif



//  if (x1 > 300 && y1 > 300 && depth > 8) return;

//  if (depth > 13) return;



  pixel_t guess[4];
/*  inv_dist_onepixel(&guess[0],encode->rec,px1,ym,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask);
  inv_dist_onepixel(&guess[1],encode->rec,px2,ym,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask);
  inv_dist_onepixel(&guess[2],encode->rec,xm,py1,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask);
  inv_dist_onepixel(&guess[3],encode->rec,xm,py2,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask); */
//  border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0);
  // order of guess: left, right, top, bottom
  

  int other=0;
  if(SHAPE_SPLITDIR(shape)) {
  // horizontal
#if HILBERT_CURVE == 1
        if(SHAPE_PIXELORDER(shape)) {
          //right-left
          if (!hasright && nn(x2,ym,px1,px2,py1,py2) && (other=1)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,1);
                output_pixel(encode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar,bdiffE,depth,0,0);
                }
          if (!hasleft  && nn2(x1,ym,px1,px2,x2,py1,py2,ym,other)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,0);
                output_pixel(encode,x1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar,bdiffE,depth,0,0);
                }
        } else {
#endif
          //left-right
          if (!hasleft  && nn(x1,ym,px1,px2,py1,py2) && (other=1)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,0);
                output_pixel(encode,x1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar,bdiffE,depth,0,0);
                }
          if (!hasright && nn2(x2,ym,px1,px2,x1,py1,py2,ym,other)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,1);
                output_pixel(encode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar,bdiffE,depth,0,0);
                }
#if HILBERT_CURVE == 1
        }
#endif
  } else {
  // vertical
#if HILBERT_CURVE == 1
        if(SHAPE_PIXELORDER(shape)) {
          //bottom-top
          if (!hasbottom && nn(xm,y2,px1,px2,py1,py2) && (other=1)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,3);
                output_pixel(encode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar,bdiffE,depth,0,1);
                }
          if (!hastop    && nn2(xm,y1,px1,px2,xm,py1,py2,y2,other)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,2);
                output_pixel(encode,xm,y1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar,bdiffE,depth,0,1);
                }
        } else {
#endif
          //top-bottom
          if (!hastop    && nn(xm,y1,px1,px2,py1,py2) && (other=1)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,2);
                output_pixel(encode,xm,y1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar,bdiffE,depth,0,1);
                }
          if (!hasbottom && nn2(xm,y2,px1,px2,xm,py1,py2,y1,other)) {
                border_guess(current_method,encode,px1,px2,py1,py2,xm,ym,guess,0,3);
                output_pixel(encode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar,bdiffE,depth,0,1);
                }
#if HILBERT_CURVE == 1
        }
#endif
  }
  depth++;
  checksplit = 1;
  
  switch(shape) {
#if HILBERT_CURVE == 0
        case SHAPE_0:
          encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft, SHAPE_0,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0,depth,checksplit); // right side
          break;
        case SHAPE_2:
          encode_recurse(encode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,   SHAPE_2,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2,depth,checksplit); // bottom side
          break;
#else
        case SHAPE_0:
          encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft, SHAPE_0A,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0B,depth,checksplit); // right side
          break;
        case SHAPE_2:
          encode_recurse(encode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,   SHAPE_2A,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2B,depth,checksplit); // bottom side
          break;
        case SHAPE_1:
          encode_recurse(encode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft,SHAPE_1A,depth,checksplit); // bottom side
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   SHAPE_1B,depth,checksplit); // top side
          break;
        case SHAPE_3:
          encode_recurse(encode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, SHAPE_3A,depth,checksplit); // right side
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,  SHAPE_3B,depth,checksplit); // left side
          break;
        case SHAPE_0A:
          encode_recurse(encode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,    SHAPE_2,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft, SHAPE_0,depth,checksplit); // bottom side
          break;
        case SHAPE_0B:
          encode_recurse(encode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft, SHAPE_0,depth,checksplit); // bottom side
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,    SHAPE_1,depth,checksplit); // top side
          break;
        case SHAPE_1A:
          encode_recurse(encode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,  SHAPE_3,depth,checksplit); // right side
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          break;
        case SHAPE_1B:
          encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,  SHAPE_0,depth,checksplit); // right side
          break;
        case SHAPE_2A:
          encode_recurse(encode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_0,depth,checksplit); // left side
          encode_recurse(encode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,  SHAPE_2,depth,checksplit); // right side
          break;
        case SHAPE_2B:
          encode_recurse(encode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,  SHAPE_2,depth,checksplit); // right side
          encode_recurse(encode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,   SHAPE_3,depth,checksplit); // left side
          break;
        case SHAPE_3A:
          encode_recurse(encode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft, SHAPE_1,depth,checksplit); // bottom side
          encode_recurse(encode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,    SHAPE_3,depth,checksplit); // top side
          break;
        case SHAPE_3B:
          encode_recurse(encode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,    SHAPE_3,depth,checksplit); // top side
          encode_recurse(encode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft, SHAPE_2,depth,checksplit); // bottom side
          break;
#endif
  }

#endif

//   for (int x = px1; x<x2; x++) pixel_set(image_pixel(encode->source,x,py1),&purple,0);
//   for (int y = py1; y<y2; y++) pixel_set(image_pixel(encode->source,px1,y),&purple,0);

//  fprintf(stderr,"done area (%i,%i)-(%i,%i)\n",x1,y1,x2,y2);
}

#if FULL_BORDERS
// recursive function for encoding
void static decode_recurse_line(coder_t *decode, int x1, int x2, int y1, int y2, int mask, int depth, uint32_t *bvar, int shape) {
  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  uint32_t p = w * h;
  int size = ilog2(p);
  int qsize_sl = quantize_log_uint32(p,QZ_SIZE_SL);
  int current_method = 0;
  assert(x2 >= x1 && y2 >= y1);
  assert(x2 == x1 || y2 == y1);
  
  if (p < 3) return;

//  fprintf(stderr,"Line from (%i,%i) to (%i,%i) at depth %i\n",x1,y1,x2,y2,depth);
  
  int method;
  int oldmask = mask;
  pixel_t *orig[4] = {
        image_pixel(decode->rec,x1,y1),               //left top
        image_pixel(decode->rec,x2,y1),               //right top
        image_pixel(decode->rec,x1,y2),               //left bottom
        image_pixel(decode->rec,x2,y2)                //right bottom
        };

  int always_split[3] = {0,0,0};
  for (int c=0; c<3; c++) {
     if (!(mask & (1 << c))) {
       if (abs(orig[0]->d[c]-orig[3]->d[c]) > ALWAYS_SPLIT_LINE_DIST*decode->qz[c]) {
                oldmask |= 1<<c;
                always_split[c] = 1;
        }                
     }
  }
        
  if (p-2 > decode->minsplitL) {
  
    mask = input_mask(decode,oldmask,qsize_sl,depth,bvar,bdiffE,1,orig[0],orig[3]);
    for (int c = 0; c < 3; c++) {
      if (!(oldmask & (1 << c))) {
        if (no_interpolation_choice(decode,x1,x2,y1,y2,~(1<<c))) {
           int i=0;
           interpolation(decode->interpolation_methods[i],decode->rec,x1,x2,y1,y2,~(1<<c),0,0,0,0,&encode->color_data);
        } else {
          if (mask & (1 << c)) {
            method = output_interpolation_choice(decode,size,c,current_method);
            interpolation(decode->interpolation_methods[method],decode->rec,x1,x2,y1,y2,~(1<<c),0,0,0,0,&encode->color_data);
            decode->method[decode->interpolation_methods[method]]++;
          }
        }
      }
    }
  }

  for (int c=0; c<3; c++) {
     if (always_split[c]) mask -= 1<<c;
  }

  if (mask == 7) return;

  pixel_t guess[4];

  int xm = (x1 + x2) / 2 ;
  int ym = (y1 + y2) / 2 ;
  assert(xm >= 0);
  assert(ym >= 0);
  
//  border_guess(current_method,decode,x1,x2,y1,y2,xm,ym,guess,0);

  depth++;
  assert(depth < DEPTH_LEVELS);
  if(x1==x2) {
  // horizontal
        border_guess(current_method,decode,x1,x2,y1,y2,xm,ym,guess,0,1);        
        input_pixel(decode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar,depth,0,0);
        if(SHAPE_PIXELORDER(shape)) {
          //bottom-top
          decode_recurse_line(decode,x1,x1,ym,y2,mask,depth,bvar,shape); 
          decode_recurse_line(decode,x1,x1,y1,ym,mask,depth,bvar,shape); 
       
        } else {
          //top-bottom
          decode_recurse_line(decode,x1,x1,y1,ym,mask,depth,bvar,shape); 
          decode_recurse_line(decode,x1,x1,ym,y2,mask,depth,bvar,shape); 
          
        }
  } else {
  // vertical
      assert(y1==y2);
      border_guess(current_method,decode,x1,x2,y1,y2,xm,ym,guess,0,3);        
      input_pixel(decode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar,depth,0,1);
        if(SHAPE_PIXELORDER(shape)) {
          //right-left
          decode_recurse_line(decode,xm,x2,y1,y1,mask,depth,bvar,shape); 
          decode_recurse_line(decode,x1,xm,y1,y1,mask,depth,bvar,shape); 
        } else {
          //left-right
          decode_recurse_line(decode,x1,xm,y1,y1,mask,depth,bvar,shape); 
          decode_recurse_line(decode,xm,x2,y1,y1,mask,depth,bvar,shape); 
        }
  }
}
#endif

// recursive function for encoding
void static decode_recurse(coder_t *decode, int x1, int x2, int y1, int y2, int mask, int previous_method, int hastop, int hasright, int hasbottom, int hasleft, shape_t shape, int depth, int checksplit) {
#if FULL_BORDERS
  int px1 = x1 - 1;
  int py1 = y1 - 1;
  int px2 = x2 + 1;
  int py2 = y2 + 1;
  int myhastop=1;
  int myhasright=1;
  int myhasbottom=1;
  int myhasleft=1;
#else
  int px1 = x1 - hasleft;
  int py1 = y1 - hastop;
  int px2 = x2 + hasright;
  int py2 = y2 + hasbottom;
  int myhastop=hastop;
  int myhasright=hasright;
  int myhasbottom=hasbottom;
  int myhasleft=hasleft;
#endif  

  if (depth == PROGRESS_BAR_VERBOSITY) fprintf(stderr,PROGRESS_BAR_CHAR);

  int w = x2 - x1 + 1;
  int h = y2 - y1 + 1;
  uint32_t p = w * h;


#ifdef DEBUGMODE
  fprintf(stderr,"%.*sarea (%i,%i)-(%i,%i) ",depth*2,"                                                                               ",x1,y1,x2,y2);
  fprintf(stderr,"shape:%s ",shapestr[shape]);
  if (hastop) fprintf(stderr,"top ");
  if (hasright) fprintf(stderr,"right ");
  if (hasbottom) fprintf(stderr,"bottom ");
  if (hasleft) fprintf(stderr,"left ");
  fprintf(stderr,"\n");
#endif

  if (px2 - px1 <= 1 && py2 - py1 <= 1) return;
  if (x1 == x2+1 || y1 == y2+1) return;

//  fprintf(stderr,"XXXX x=%i y=%i  px=%i py=%i  depth=%i  shape=%s \n",x2-x1+1,y2-y1+1,px2-px1+1,py2-py1+1,depth,shapestr[shape]);

  assert(x2 >= x1 && y2 >= y1);



//  int method;
  int current_method = previous_method;

  uint32_t bvar[3] = {};
  uint32_t bdiffE[3] = {};
//  pixel_t firstguess;
  compute_bvar(decode,x1,y1,x2,y2,0,&bvar[0],myhastop,myhasright,myhasbottom,myhasleft); //,&firstguess);
  compute_bdiffE(decode,x1,y1,x2,y2,0,&bdiffE[0],myhastop,myhasright,myhasbottom,myhasleft); 

  pixel_t *orig[4] = {
        image_pixel(decode->rec,px1,py1),               //left top
        image_pixel(decode->rec,px2,py1),               //right top
        image_pixel(decode->rec,px1,py2),               //left bottom
        image_pixel(decode->rec,px2,py2)                //right bottom
        };



  
  int oldmask = mask;

#if ALWAYS_SPLIT_HEURISTIC
  int always_split[3] = {0,0,0};
  for (int c=0; c<3; c++) {
     if (!(mask & (1 << c))) {
       if (bvar[c] > ALWAYS_SPLIT_AREA_DIST*decode->qz[c]) {
                oldmask |= 1<<c;
                always_split[c] = 1;
        }                
     }
  }
#endif

  
  if (p > decode->minsplitA && depth%SPLIT_MODULO == 0) {

//           interpolation(decode->interpolation_methods[0],decode->rec,x1,x2,y1,y2,oldmask,myhastop,myhasright,myhasbottom,myhasleft,&decode->color_data);
  
  int qsize = quantize_log_uint32(p,QZ_SIZE);
//  int qsize_s = quantize_log_uint32(p,QZ_SIZE_S+1)-1;
    mask = input_mask(decode,oldmask,qsize,depth,bvar,bdiffE,0,NULL,NULL);
/*
    for (int c = 0; c < 3; c++) {
      if (!(oldmask & (1 << c))) {
        if (no_interpolation_choice(decode,x1,x2,y1,y2,~(1<<c))) {
           int i=0;
           interpolation(decode->interpolation_methods[i],decode->rec,x1,x2,y1,y2,~(1<<c),myhastop,myhasright,myhasbottom,myhasleft,&decode->color_data);
        } else {
          if (mask & (1 << c)) {
            method = input_interpolation_choice(decode,qsize,c,current_method);
            interpolation(decode->interpolation_methods[method],decode->rec,x1,x2,y1,y2,~(1<<c),myhastop,myhasright,myhasbottom,myhasleft,&decode->color_data);
            decode->method[decode->interpolation_methods[method]]++;
          }
        }
      }
    }
*/    
  }

#if ALWAYS_SPLIT_HEURISTIC
  for (int c=0; c<3; c++) {
     if (always_split[c]) mask -= 1<<c;
  }
#endif

  if (mask != oldmask)
      interpolation(decode->interpolation_methods[0],decode->rec,x1,x2,y1,y2,oldmask,myhastop,myhasright,myhasbottom,myhasleft,&decode->color_data);

  if (mask == 7) return;

  int xm = (x1 + x2) / 2 ;
  int ym = (y1 + y2) / 2 ;
#if HILBERT_CURVE == 0
  shape = decide_splitdir(decode->rec,x1,x2,y1,y2,px1,px2,py1,py2,mask,&xm,&ym); 
/*
  if (x2-x1 > y2-y1) {          //                      |
        shape = SHAPE_0;        // vertical split     L | R
                                //                      |
                                
  } else {                      //                      T
        shape = SHAPE_2;        // horizontal split   -----
  }                             //                      B
*/
#endif 


#if FULL_BORDERS

  if ( x1==x2 && y1==y2) {
    // one pixel

    wpixel_t myguess = {};
    wpixel_add(&myguess,image_pixel(decode->rec,x1,py1),2,0);
    wpixel_add(&myguess,image_pixel(decode->rec,x1,py2),2,0);
    wpixel_add(&myguess,image_pixel(decode->rec,px1,y1),2,0);
    wpixel_add(&myguess,image_pixel(decode->rec,px2,y1),2,0);
    wpixel_add(&myguess,image_pixel(decode->rec,px1,py1),-1,0);
    wpixel_add(&myguess,image_pixel(decode->rec,px2,py1),-1,0);
    wpixel_add(&myguess,image_pixel(decode->rec,px1,py2),-1,0);
    wpixel_add(&myguess,image_pixel(decode->rec,px2,py2),-1,0);
    wpixel_div(&myguess,4,0);
    pixel_t grad;
    wpixel_set(&grad,&myguess,0);
    pixel_t *pixels[5] = {
            image_pixel(decode->rec,x1,py1),
            image_pixel(decode->rec,x1,py2),
            image_pixel(decode->rec,px1,y1),
            image_pixel(decode->rec,px2,y1),
            &grad };
    pixel_t guess;
//    int weights[5] = {1,1,1,1,1};
//    pixel_med(&guess, pixels, weights, 5, mask);
    for (int c=0; c<3; c++) {
        guess.d[c] = median5(pixels[0]->d[c],pixels[1]->d[c],pixels[2]->d[c],pixels[3]->d[c],pixels[4]->d[c]);
    }
    input_pixel(decode,x1,y1,&guess,p,mask,orig[0],orig[1],orig[2],orig[3],bvar,bdiffE,depth,1,0);
    return;
  }

  int use_firstguess = 1;
  if (w < h-2 || w > h+2)   use_firstguess = 0;  
  // if (w == 1 && h == 1) use_firstguess = 1;

  if (use_firstguess && USE_BORDER_GUESSING) {
  input_pixel(decode,xm,ym,&firstguess,p,mask,orig[2],orig[1],orig[0],orig[3],bvar,depth,1,0);
  if(SHAPE_SPLITDIR(shape)) {
  // horizontal
        if(SHAPE_PIXELORDER(shape)) {
          //right-left
          decode_recurse_line(decode,xm,px2,ym,ym,mask,depth,bvar,shape); 
          decode_recurse_line(decode,px1,xm,ym,ym,mask,depth,bvar,shape); 
        } else {
          //left-right
          decode_recurse_line(decode,px1,xm,ym,ym,mask,depth,bvar,shape); 
          decode_recurse_line(decode,xm,px2,ym,ym,mask,depth,bvar,shape); 
        }
  } else {
  // vertical
        if(SHAPE_PIXELORDER(shape)) {
          //bottom-top
          decode_recurse_line(decode,xm,xm,ym,py2,mask,depth,bvar,shape); 
          decode_recurse_line(decode,xm,xm,py1,ym,mask,depth,bvar,shape); 
        } else {
          //top-bottom
          decode_recurse_line(decode,xm,xm,py1,ym,mask,depth,bvar,shape); 
          decode_recurse_line(decode,xm,xm,ym,py2,mask,depth,bvar,shape); 
        }
  }
  
  } else {
    if(SHAPE_SPLITDIR(shape)) {
    // horizontal
          decode_recurse_line(decode,px1,px2,ym,ym,mask,depth,bvar,shape);
    } else {
    // vertical
          decode_recurse_line(decode,xm,xm,py1,py2,mask,depth,bvar,shape);
    }
  }
  depth++;
  checksplit = 1;
  
  switch(shape) {
        case SHAPE_0:
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_0A,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0B,depth,checksplit); // right side
          break;
        case SHAPE_1:
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,0,hasright,hasbottom,hasleft,  SHAPE_1A,depth,checksplit); // bottom side
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   SHAPE_1B,depth,checksplit); // top side
          break;
        case SHAPE_2:
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,0,hasleft,     SHAPE_2A,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2B,depth,checksplit); // bottom side
          break;
        case SHAPE_3:
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, SHAPE_3A,depth,checksplit); // right side
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,SHAPE_3B,depth,checksplit); // left side
          break;
        case SHAPE_0A:
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,0,hasleft,     SHAPE_2,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_0,depth,checksplit); // bottom side
          break;
        case SHAPE_0B:
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,0,hasright,hasbottom,hasleft,  SHAPE_0,depth,checksplit); // bottom side
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   SHAPE_1,depth,checksplit); // top side
          break;
        case SHAPE_1A:
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, SHAPE_3,depth,checksplit); // right side
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,SHAPE_1,depth,checksplit); // left side
          break;
        case SHAPE_1B:
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0,depth,checksplit); // right side
          break;
        case SHAPE_2A:
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_0,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_2,depth,checksplit); // right side
          break;
        case SHAPE_2B:
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, SHAPE_2,depth,checksplit); // right side
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,SHAPE_3,depth,checksplit); // left side
          break;
        case SHAPE_3A:
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,0,hasright,hasbottom,hasleft,  SHAPE_1,depth,checksplit); // bottom side
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   SHAPE_3,depth,checksplit); // top side
          break;
        case SHAPE_3B:
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,0,hasleft,     SHAPE_3,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2,depth,checksplit); // bottom side
          break;
  }
#else


#if HILBERT_CURVE == 1
  xm = (px1+px2)/2;
  ym = (py1+py2)/2;
#endif

  pixel_t guess[4];

/*  inv_dist_onepixel(&guess[0],decode->rec,px1,ym,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask);
  inv_dist_onepixel(&guess[1],decode->rec,px2,ym,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask);
  inv_dist_onepixel(&guess[2],decode->rec,xm,py1,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask);
  inv_dist_onepixel(&guess[3],decode->rec,xm,py2,px1,px2,py1,py2,x1,x2,y1,y2,px1,px2,py1,py2,0*mask); */

//  border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0);
  int other=0;
  // order of guess: left, right, top, bottom
  if(SHAPE_SPLITDIR(shape)) {
  // horizontal
#if HILBERT_CURVE == 1
        if(SHAPE_PIXELORDER(shape)) {
          //right-left
          if (!hasright && nn(x2,ym,px1,px2,py1,py2) && (other=1)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,1);input_pixel(decode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar,bdiffE,depth,0,0);}
          if (!hasleft  && nn2(x1,ym,px1,px2,x2,py1,py2,ym,other)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,0);input_pixel(decode,x1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar,bdiffE,depth,0,0);}
        } else {
#endif
          //left-right
          if (!hasleft  && nn(x1,ym,px1,px2,py1,py2) && (other=1)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,0);input_pixel(decode,x1,ym,&guess[0],p,mask,orig[1],orig[3],orig[0],orig[2],bvar,bdiffE,depth,0,0);}
          if (!hasright && nn2(x2,ym,px1,px2,x1,py1,py2,ym,other)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,1);input_pixel(decode,x2,ym,&guess[1],p,mask,orig[2],orig[0],orig[3],orig[1],bvar,bdiffE,depth,0,0);}
#if HILBERT_CURVE == 1
        }
#endif
  } else {
  // vertical
#if HILBERT_CURVE == 1
        if(SHAPE_PIXELORDER(shape)) {
          //bottom-top
          if (!hasbottom && nn(xm,y2,px1,px2,py1,py2) && (other=1)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,3);input_pixel(decode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar,bdiffE,depth,0,1);}
          if (!hastop    && nn2(xm,y1,px1,px2,xm,py1,py2,y2,other)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,2);input_pixel(decode,xm,y1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar,bdiffE,depth,0,1);}
        } else {
#endif
          //top-bottom
          if (!hastop    && nn(xm,y1,px1,px2,py1,py2) && (other=1)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,2);input_pixel(decode,xm,y1,&guess[2],p,mask,orig[3],orig[2],orig[1],orig[0],bvar,bdiffE,depth,0,1);}
          if (!hasbottom && nn2(xm,y2,px1,px2,xm,py1,py2,y1,other)) {border_guess(current_method,decode,px1,px2,py1,py2,xm,ym,guess,0,3);input_pixel(decode,xm,y2,&guess[3],p,mask,orig[0],orig[1],orig[2],orig[3],bvar,bdiffE,depth,0,1);}
#if HILBERT_CURVE == 1
        }
#endif
  }
  depth++;
  checksplit = 1;
  
  switch(shape) {
#if HILBERT_CURVE == 0
        case SHAPE_0:
          decode_recurse(decode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft, SHAPE_0,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0,depth,checksplit); // right side
          break;
        case SHAPE_2:
          decode_recurse(decode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,   SHAPE_2,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2,depth,checksplit); // bottom side
          break;
#else
        case SHAPE_0:
          decode_recurse(decode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft, SHAPE_0A,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,SHAPE_0B,depth,checksplit); // right side
          break;
        case SHAPE_2:
          decode_recurse(decode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,   SHAPE_2A,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft,SHAPE_2B,depth,checksplit); // bottom side
          break;
        case SHAPE_1:
          decode_recurse(decode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft,SHAPE_1A,depth,checksplit); // bottom side
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,   SHAPE_1B,depth,checksplit); // top side
          break;
        case SHAPE_3:
          decode_recurse(decode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0, SHAPE_3A,depth,checksplit); // right side
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,  SHAPE_3B,depth,checksplit); // left side
          break;
        case SHAPE_0A:
          decode_recurse(decode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,    SHAPE_2,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft, SHAPE_0,depth,checksplit); // bottom side
          break;
        case SHAPE_0B:
          decode_recurse(decode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft, SHAPE_0,depth,checksplit); // bottom side
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,    SHAPE_1,depth,checksplit); // top side
          break;
        case SHAPE_1A:
          decode_recurse(decode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,  SHAPE_3,depth,checksplit); // right side
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          break;
        case SHAPE_1B:
          decode_recurse(decode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_1,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,  SHAPE_0,depth,checksplit); // right side
          break;
        case SHAPE_2A:
          decode_recurse(decode,x1,xm,y1,y2,mask,current_method,hastop,0,hasbottom,hasleft,   SHAPE_0,depth,checksplit); // left side
          decode_recurse(decode,xm+1,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,1,  SHAPE_2,depth,checksplit); // right side
          break;
        case SHAPE_2B:
          decode_recurse(decode,xm,x2,y1,y2,mask,current_method,hastop,hasright,hasbottom,0,  SHAPE_2,depth,checksplit); // right side
          decode_recurse(decode,x1,xm-1,y1,y2,mask,current_method,hastop,1,hasbottom,hasleft,   SHAPE_3,depth,checksplit); // left side
          break;
        case SHAPE_3A:
          decode_recurse(decode,x1,x2,ym,y2,mask,current_method,0,hasright,hasbottom,hasleft, SHAPE_1,depth,checksplit); // bottom side
          decode_recurse(decode,x1,x2,y1,ym-1,mask,current_method,hastop,hasright,1,hasleft,    SHAPE_3,depth,checksplit); // top side
          break;
        case SHAPE_3B:
          decode_recurse(decode,x1,x2,y1,ym,mask,current_method,hastop,hasright,0,hasleft,    SHAPE_3,depth,checksplit); // top side
          decode_recurse(decode,x1,x2,ym+1,y2,mask,current_method,1,hasright,hasbottom,hasleft, SHAPE_2,depth,checksplit); // bottom side
          break;
#endif
  }

#endif


//  fprintf(stderr,"done area (%i,%i)-(%i,%i)\n",x1,y1,x2,y2);
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
        if (nb_methods > 1) fprintf(stderr,"Nb of methods: %i \n",enc->nb_interpolation_methods);
}

void write_context_trees(coder_t *enc) {
    for(int c=0; c<3; c++) {
       if (enc->phase == 2) {
         output_context_tree_p2(enc->ctx->diff_p2[c],enc,0,c);
       } else {
         output_context_tree(enc->ctx->diff[c],enc,0,c);
       }
       fprintf(stderr,"Wrote context tree %i  (%g bytes)\n", c, enc->ctx->bitsContextTreeData[c]*log4k_base/8);
    }
}
void read_context_trees(coder_t *enc) {
    for(int c=0; c<3; c++) {
       fprintf(stderr,"%i ", c);
       input_context_tree_p2(enc->ctx->diff_p2[c],enc,0,c);
    }
    fprintf(stderr,"\n");
}

//#define CHART_WIDTH 8192
//#define CHART_HEIGHT 8192
//#define CHART_SEP 2

#define CHART_WIDTH 1024
#define CHART_HEIGHT 1024
#define CHART_SEP 0

void actual_make_auto_indexing_chart(int only_used_colors, colors *colors) {
  image_t rec = {};
  image_init(&rec,CHART_WIDTH,CHART_HEIGHT,&purple);
  for(int y=0; y<CHART_HEIGHT; y++) {
    for(int x=0; x<CHART_WIDTH; x++) {

      int bucket_i = (16*y*510/(CHART_HEIGHT-1))/(16*SIZE_I);
      int bucket_q = (16*x*510/(CHART_WIDTH-1))/(16*SIZE_Q);
      int bucket_ip = (16*y*510/(CHART_HEIGHT-1))%(16*SIZE_I);
      int bucket_qp = (16*x*510/(CHART_WIDTH-1))%(16*SIZE_Q);
//      int bucket_y = (-8*SIZE_Y + bucket_qp*BUCKETS_Q + bucket_ip*BUCKETS_I)/(64*SIZE_Y);
//      int bucket_yp = ((-8*SIZE_Y + bucket_qp*BUCKETS_Q + bucket_ip*BUCKETS_I))%(64*SIZE_Y);
      int bucket_y = (bucket_ip*BUCKETS_I)/(32*SIZE_Y);
      int bucket_yp = (bucket_ip*BUCKETS_I)%(32*SIZE_Y);

                //(bucket_qp*BUCKETS_Q + bucket_ip*BUCKETS_I)/(4*SIZE_Y);
//      int bucket_y = (bucket_qp*BUCKETS_Q/2/SIZE_Y + bucket_ip*BUCKETS_I/2/SIZE_Y)/2;
//      int bucket_y = (((x*510/(CHART_WIDTH-1))%SIZE_Q)*BUCKETS_Q/2  ) / SIZE_Y;
        //1;//(bucket_qp*BUCKETS_Q + bucket_ip*BUCKETS_I)/4  %SIZE_Y;
        //1; //((x*510/(CHART_WIDTH-1))%SIZE_Q)*BUCKETS_Q/2 % SIZE_Y;

/*      int bucket_y = (y*255/(CHART_HEIGHT-1))/SIZE_Y;
      int bucket_i = (x*510/(CHART_WIDTH-1))/SIZE_I;
      int bucket_q = ((x*510/(CHART_WIDTH-1))%SIZE_I)*BUCKETS_I / SIZE_Q;
      int bucket_yp = (y*255/(CHART_HEIGHT-1))%SIZE_Y ;
      int bucket_ip = (x*510/(CHART_WIDTH-1))%SIZE_I ;
      int bucket_qp = ((x*510/(CHART_WIDTH-1))%SIZE_I)*BUCKETS_I % SIZE_Q  ;
  */
      pixel_t * p = image_pixel(&rec,x,y);
      if (bucket_exists(bucket_y,bucket_i,bucket_q)
            && (only_used_colors==0 || colors->bucketYIQ[bucket_y][bucket_i][bucket_q].count)
            && bucket_yp>CHART_SEP*64-1 && bucket_ip>CHART_SEP*8-1 && bucket_qp>CHART_SEP*8-1) {
        if (only_used_colors && colors->bucketYIQ[bucket_y][bucket_i][bucket_q].count < MAX_PER_BUCKET) {
          int i = (bucket_qp*MAX_PER_BUCKET)/(16*SIZE_Q);
          if (i>=0 && i < colors->bucketYIQ[bucket_y][bucket_i][bucket_q].count) {
            p->d[0] = STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[i].d[0]);
            p->d[1] = STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[i].d[1]);
            p->d[2] = STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[i].d[2]);
          } else {
            p->d[0] = STRETCH(0);
            p->d[1] = STRETCH(255);
            p->d[2] = STRETCH(255);
          }
        } else {
          p->d[0] = STRETCH((bucket_y)*SIZE_Y+((bucket_yp/32)%SIZE_Y));
          p->d[1] = STRETCH((bucket_i)*SIZE_I+((bucket_qp/16)%SIZE_I));
          p->d[2] = STRETCH((bucket_q)*SIZE_Q+((bucket_qp/16)%SIZE_Q));
//        p->d[0] = STRETCH((2*bucket_y+1)*SIZE_Y/2);
//        p->d[1] = STRETCH((2*bucket_i+1)*SIZE_I/2);
//        p->d[2] = STRETCH((2*bucket_q+1)*SIZE_Q/2);
          if (only_used_colors) {
                //round_color(colors,p);
              if (p->d[0] < STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[0].d[0])
              ||  p->d[1] < STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[0].d[1])
              ||  p->d[2] < STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[0].d[2])
              ||  p->d[0] > STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[1].d[0])
              ||  p->d[1] > STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[1].d[1])
              ||  p->d[2] > STRETCH(colors->bucketYIQ[bucket_y][bucket_i][bucket_q].color[1].d[2])) {
                p->d[0] = STRETCH(255);
                p->d[1] = STRETCH(255);
                p->d[2] = STRETCH(255);
              }
                
          }
        }
      } else {
        p->d[0] = STRETCH(0);
        p->d[1] = STRETCH(255);
        p->d[2] = STRETCH(255);
      }
    }
  }
  yiq2rgb(&rec);
  image_save("auto_indexing_chart.pnm",&rec);
  image_free(&rec);
  fprintf(stderr,"Saved auto-indexing chart to auto_indexing_chart.pnm\n");
}

void make_auto_indexing_chart() {
    actual_make_auto_indexing_chart(0,NULL);
}

void make_auto_indexing_chart_used(colors *colors) {
    actual_make_auto_indexing_chart(1,colors);
}

void progress_bar_firstline() {
    fprintf(stderr," 0%%");
    for (int i=1; i < (1<<PROGRESS_BAR_VERBOSITY)-3 ; i++) {
        if (i == (1<<(PROGRESS_BAR_VERBOSITY-1))-1) {
               fprintf(stderr,"50%%");
               i += 2;
        } else if (i == (1<<(PROGRESS_BAR_VERBOSITY-2))-1) {
               fprintf(stderr,"25%%");
               i += 2;
        } else if (i == (3<<(PROGRESS_BAR_VERBOSITY-2))-1) {
               fprintf(stderr,"75%%");
               i += 2;
        } else if (i % (1<<(PROGRESS_BAR_VERBOSITY-3)) == 0) {
               fprintf(stderr,".");
        } else {
               fprintf(stderr," ");
        }
    }
    fprintf(stderr,"100%%\n [");
}


// encode an image
void static encode(image_t *img, double epsilon, double epsiloni, double epsilonq, double factor, double factor_s, double qzy, double qzi, double qzq, char *out, char *methods, int skipsplitA, int skipsplitL, int cutoff, int *sb, int t_m, int phase, int auto_indexing_chart) {
  image_t rec = {};
  image_init(&rec,img->w,img->h,&purple);
  if (qzy < 1.0) qzy = 1.0;
  if (qzi < 1.0) qzi = 1.0;
  if (qzq < 1.0) qzq = 1.0;
  if (qzy > 255.0) qzy = 255.0;
  if (qzi > 255.0) qzi = 255.0;
  if (qzq > 255.0) qzq = 255.0;
  int msA=skipsplitA;
  int msL=skipsplitL;
  ctx_t *ctx=malloc(sizeof(ctx_t));
  coder_t enc = {.source = img,.rec = &rec,.maxD = {epsilon,epsiloni,epsilonq},.factor = factor, .factor_s = factor_s, .traversal_method = t_m,
                 .outputted_pixels = {0,0,0}, .qz={(int)(qzy*COLOR_STRETCH+0.5), (int)(qzi*COLOR_STRETCH+0.5), (int)(qzq*COLOR_STRETCH+0.5)}, .pLeft=0, .pRight=0, .pTop=0, .pBottom=0,
                 .sb = {sb[0], sb[1], sb[2]},
                 .avg_Lsplit_depth = {},
                 .avg_Asplit_depth = {}, .phase = phase,
                 .color_data = {}, .nb_pixels = img->w * img->h,
                 .method={}, .method_gain={}, .minsplitA=msA, .minsplitL=msL, .ctx=ctx, .range_min={}, .range_max={}};
  parse_interpolation_methods(&enc, methods);
  fprintf(stderr,"Encoding %ix%i image, ",(int)(img->w),(int)(img->h));
  if (phase>0) fprintf(stderr,"phase %i, ", phase);
  if (t_m == 0) fprintf(stderr,"mode: splitting\n");
  if (t_m == 1) fprintf(stderr,"mode: FFV1\n");

  if (enc.maxD[0] == 0 && enc.qz[0] == COLOR_STRETCH
  && enc.maxD[1] == 0 && enc.qz[1] == COLOR_STRETCH
  && enc.maxD[2] == 0 && enc.qz[2] == COLOR_STRETCH ) {
        fprintf(stderr,"Lossless compression on all channels.\n");
  } else {
        fprintf(stderr,"Lossy compression settings:\n");
    for(int c=0; c<3; c++) {
        fprintf(stderr,"  %c: ",planesYIQ[c]);
        if (enc.maxD[c] == 0 && enc.qz[c] == COLOR_STRETCH) {
           fprintf(stderr,"lossless");
        } else {
           if (t_m == 0) {
                fprintf(stderr,"pixel dist < %i or area dist < %.1f and avg dist < %.1f, qz=%i",enc.maxD[c],enc.maxD[c]*factor,enc.maxD[c]*factor_s,enc.qz[c]/COLOR_STRETCH);
           } else {
                fprintf(stderr,"quantization=%i",enc.qz[c]/COLOR_STRETCH);
           }
        }
        fprintf(stderr,"\n");
//        if (enc.qz[c]>COLOR_STRETCH && enc.maxD[c]*COLOR_STRETCH < 2*enc.qz[c])
        if (enc.qz[c]>COLOR_STRETCH && enc.maxD[c]*COLOR_STRETCH < enc.qz[c])
                fprintf(stderr,"    Warning: quantization on channel %c is %i while max distance is only %i!\n",planesYIQ[c],enc.qz[c]/COLOR_STRETCH,enc.maxD[c]);
    }
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

  if (cutoff != 8)
  fprintf(stderr,"Using RAC range [%g-%g] (cutoff %i)\n",cutoff/4096.0,(4096-cutoff)/4096.0,cutoff);

  if (t_m == 0) {
    if (msA>0) {
        if (msA < img->w * img->h) {
                fprintf(stderr,"Areas smaller than %i pixels will split immediately.\n",msA);
        } else {
                fprintf(stderr,"All areas will split immediately.\n");
        }
    } else {
        fprintf(stderr,"No areas will split immediately.\n");
    }
  }
  if (msL>0 && t_m == 0) fprintf(stderr,"Lines shorter than %i pixels will always be split\n",msL);
  fputc((skipsplitA >> 0) & 0xFF,f);
  fputc((skipsplitA >> 8) & 0xFF,f);
  fputc((skipsplitA >> 16) & 0xFF,f);
  fputc((skipsplitA >> 24) & 0xFF,f);
  fputc((skipsplitL >> 0) & 0xFF,f);
  fputc((skipsplitL >> 8) & 0xFF,f);

  fputc((sb[0] >> 0) & 0xFF,f);
  fputc((sb[1] >> 0) & 0xFF,f);
  fputc((sb[2] >> 0) & 0xFF,f);
  if (sb[0] < 8 || sb[1] < 9 || sb[2] < 9) 
          fprintf(stderr,"Significant mantissa bits: Y: %i/8  I: %i/9  Q: %i/9\n", sb[0], sb[1], sb[2]);


  ctx_init(enc.ctx,cutoff,phase);

  if (phase == 2) {
    fprintf(stderr,"Trying to read context_trees.dat\n");
    FILE *g = fopen("context_trees.dat","rb");
    symb_init_read(&enc.coder,g,cutoff);
    fprintf(stderr,"Reading context trees: ");
    read_context_trees(&enc);
    fclose(g);
    ctx_init(enc.ctx,cutoff,phase);
  }


  symb_init_write(&enc.coder,f,cutoff);
//  extra_ctx_init(&enc,cutoff);



  output_interpolation_methods(&enc);
  for (int i = 0; i < 16; i++) {
    symb_put_simple_bit(&enc.coder,(imagic >> i) & 1,NULL);
  }

  image_set_alpha(img,0);
  rgb2yiq(img);

  int use_color_data = 1;
  if (enc.qz[0] > COLOR_STRETCH) use_color_data = 0;   //*SIZE_Y/3
  if (enc.qz[1] > COLOR_STRETCH) use_color_data = 0; //*SIZE_I/3
  if (enc.qz[2] > COLOR_STRETCH) use_color_data = 0;   //*SIZE_Q/3
   
  if (auto_indexing_chart == 1) { fprintf(stderr,"Doing auto-indexing. "); use_color_data = 1; }
  else if (auto_indexing_chart == 2) { fprintf(stderr,"Not doing auto-indexing. "); use_color_data = 0; }
  else if (use_color_data) fprintf(stderr,"Considering auto-indexing. ");
  else fprintf(stderr,"Not considering auto-indexing. ");
  
  fprintf(stderr,"Scanning colors used in image...\n");
  int min[3] = {255,510,510};
  int max[3] = {0,0,0};
  for(int y = 0; y < img->h; y++) {
     for(int x = 0; x < img->w; x++) {
        pixel_t *p = image_pixel(img,x,y);
        if (use_color_data) add_color(&enc.color_data, p);
        for(int c=0; c<3; c++) {
            if (UNSTRETCH(p->d[c]) < min[c]) min[c]=UNSTRETCH(p->d[c]);
            if (UNSTRETCH(p->d[c]) > max[c]) max[c]=UNSTRETCH(p->d[c]);
        }
     }
  }
  fprintf(stderr,"Color ranges: ");
  for(int c=0; c<3; c++) {
        fprintf(stderr,"%c in [%i,%i]   ",planesYIQ[c],min[c],max[c]);
        enc.range_min[c]=min[c];
        enc.range_max[c]=max[c];
        output_min_max(&enc,enc.range_min[c],enc.range_max[c],0,(c==0?255:510),c);
  }
  fprintf(stderr,"\n");
  
  if (use_color_data && !auto_indexing_chart) {
    if (!auto_indexing_heuristic(&enc.color_data)) {
        fprintf(stderr,"reducing color table.\n");
        convert_indexed_clusters_to_full_buckets(&enc.color_data);
        if (!auto_indexing_heuristic(&enc.color_data)) {
           use_color_data=0;
        }
        if (use_color_data == 0) fprintf(stderr,"not using auto-indexing.\n");
    }
//  } else {
//    if (!use_color_data && auto_indexing_chart == 0) fprintf(stderr,"Coarse quantization, so ");
  } else {
    if (use_color_data && epsilon > 10) {
        fprintf(stderr,"Very lossy compression, reducing color table.\n");
        convert_indexed_clusters_to_full_buckets(&enc.color_data);
    }
  }
  if (use_color_data) {
    enc.color_data.used = 1;
    if (!auto_indexing_chart) fprintf(stderr,"using auto-indexing. ");
    fprintf(stderr,"Outputting color table.\n");
    symb_put_simple_bit(&enc.coder,1,NULL);
    output_colors(&enc);
//    fprintf(stderr,"Using auto-indexing. Wrote color table, took %g bytes.\n", enc.ctx->bitsColorData*log4k_base/8);
    fprintf(stderr,"Compressed color table size: %i bytes\n",
        (int) ((enc.ctx->bitsColorData[0]+enc.ctx->bitsColorData[1]+enc.ctx->bitsColorData[2]+enc.ctx->bitsColorData_count+enc.ctx->bitsColorData_Y)   *log4k_base/8 +0.5));
    fprintf(stderr,"   (Y-vector:%i  YIQ-counts:%i  Y:%i  I:%i  Q:%i)\n", 
        (int) (enc.ctx->bitsColorData_Y*log4k_base/8+0.5),
        (int) (enc.ctx->bitsColorData_count*log4k_base/8+0.5),
        (int) (enc.ctx->bitsColorData[0]*log4k_base/8+0.5),
        (int) (enc.ctx->bitsColorData[1]*log4k_base/8+0.5),
        (int) (enc.ctx->bitsColorData[2]*log4k_base/8+0.5)        );
    enc.color_data.max_diff[0] = enc.qz[0]/COLOR_STRETCH ;
    enc.color_data.max_diff[1] = enc.qz[1]/COLOR_STRETCH ;
    enc.color_data.max_diff[2] = enc.qz[2]/COLOR_STRETCH ;
//    fprintf(stderr,"Max color diffs: %i %i %i\n", enc.color_data.max_diff[0], enc.color_data.max_diff[1], enc.color_data.max_diff[2]);
  } else {
//    if (auto_indexing_chart == 0) fprintf(stderr,"not using auto-indexing.\n");
    symb_put_simple_bit(&enc.coder,0,NULL);
  }
  if (auto_indexing_chart == 1) make_auto_indexing_chart_used(&enc.color_data);

  int y1 = 0, y2 = img->h-1, x1 = 0, x2 = img->w-1;


  if (phase == 2) {
    symb_put_simple_bit(&enc.coder,1,NULL);
    write_context_trees(&enc);
  } else {
    symb_put_simple_bit(&enc.coder,0,NULL);
  }

  fprintf(stderr,"\nHeader written, starting actual encoding.\n");
  progress_bar_firstline();
  if (t_m == 1) {
    symb_put_simple_bit(&enc.coder,1,NULL);
    encode_ffv1(&enc,x1,x2,y1,y2);
  } else {
    symb_put_simple_bit(&enc.coder,0,NULL);
    uint32_t p= img->w * img->h;
    uint32_t bvar[3]={};
    output_pixel(&enc,x1,y1,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
    output_pixel(&enc,x2,y1,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
    output_pixel(&enc,x2,y2,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
    output_pixel(&enc,x1,y2,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
#if FULL_BORDERS
    int mask=0;
    encode_recurse_line(&enc,x1,x2,y1,y1,mask,0,bvar,bvar,0);
    encode_recurse_line(&enc,x1,x2,y2,y2,mask,0,bvar,bvar,0);
    encode_recurse_line(&enc,x1,x1,y1,y2,mask,0,bvar,bvar,0);
    encode_recurse_line(&enc,x2,x2,y1,y2,mask,0,bvar,bvar,0);
    encode_recurse(&enc,x1+1,x2-1,y1+1,y2-1,0,0,0,0,0,0,SHAPE_0,0,1);
#else
    encode_recurse(&enc,x1,x2,y1,y2,0,0,0,0,0,0,SHAPE_0,0,1);
#endif
  }

  for (int i = 0; i < 16; i++) {
    symb_put_simple_bit(&enc.coder,(imagic >> i) & 1,NULL);
  }
  uint32_t crc=image_crc(&rec);
  for (int i = 0; i < 32; i++) {
    symb_put_simple_bit(&enc.coder,(crc >> i) & 1,NULL);
  }

  symb_flush(&enc.coder);
  fclose(f);

  fprintf(stderr,"]\n\n");


  yiq2rgb(&rec);
  image_save("debug_output_rec.pnm",&rec);
//  yiq2rgb(img);
//  image_save("debug_output_src.pnm",img);

  if (phase < 2) {
  fprintf(stderr,"Context tree sizes:\n");
  for(int c = 0; c<3 ; c++) {
   int contexts_split=0;
   int contexts_used=0;
   for (int i=0; i<=enc.ctx->treesize[0][c] ; i++) {
     if (enc.ctx->diff[c][i].property == 0) {
         contexts_used++;
     } else {
         contexts_split++;
     }
   }
  if (contexts_split>0)   fprintf(stderr,"  %c pixel context tree:  %i nodes, %i leafs, %lu kB\n",planesYIQ[c],contexts_split,contexts_used,(unsigned long) contexts_used * sizeof(context_leaf)/1024);
  }
  }
#if FULL_BORDERS
  for(int c = 0; c<3 ; c++) {
   int contexts_split=0;
   int contexts_used=0;
   for (int i=0; i<=enc.ctx->treesize[3][c] ; i++) {
     if (enc.ctx->diff1[c][i].property == 0) {
         contexts_used++;
     } else {
         contexts_split++;
     }
   }
  if (contexts_split>0)   fprintf(stderr,"Last pixel context tree (%c) nodes: %i, leafs: %i (%lu bytes)\n",planesYIQ[c],contexts_split,contexts_used,(unsigned long) contexts_used * sizeof(context_leaf));
  }
  for(int c = 0; c<3 ; c++) {
   int contexts_split=0;
   int contexts_used=0;
   for (int i=0; i<=enc.ctx->treesize[1][c] ; i++) {
     if (enc.ctx->splitContL[c][i].property == 0) {
         contexts_used++;
     } else {
         contexts_split++;
     }
   }
  if (contexts_split>0) fprintf(stderr,"Line split context tree (%c) nodes: %i, leafs: %i (%lu bytes)\n",planesYIQ[c],contexts_split,contexts_used,(unsigned long) contexts_used * sizeof(context_leaf_bit));
  }
#endif  
  for(int c = 0; c<3 ; c++) {
   int contexts_split=0;
   int contexts_used=0;
   for (int i=0; i<=enc.ctx->treesize[2][c] ; i++) {
     if (enc.ctx->splitContA[c][i].property == 0) {
         contexts_used++;
     } else {
         contexts_split++;
     }
   }
  if (contexts_split>0)   fprintf(stderr,"  %c area split context tree: %i nodes, %i leafs, %lu kB\n",planesYIQ[c],contexts_split,contexts_used,(unsigned long) contexts_used * sizeof(context_leaf_bit)/1024);
  }
  if (enc.splits[0][0]+enc.splits[1][0]+enc.splits[2][0] > 0) {
          fprintf(stderr,"Area Splits: Y:%i  I:%i  Q:%i\n",enc.splits[0][0],enc.splits[1][0],enc.splits[2][0]);
  }
  fprintf(stderr,"Asymmetric splits: %i\n",asym_splits);
  fprintf(stderr,"Outputted pixels: %i(%.2f%%) %i(%.2f%%) %i(%.2f%%)\n",enc.outputted_pixels[0],100.0 * enc.outputted_pixels[0] / (img->w * img->h),enc.outputted_pixels[1],100.0 * enc.outputted_pixels[1] / (img->w * img->h),enc.outputted_pixels[2],100.0 * enc.outputted_pixels[2] / (img->w * img->h));
#if FULL_BORDERS
  fprintf(stderr,"Line Splits: Y:%i  I:%i  Q:%i\n",enc.splits[0][1],enc.splits[1][1],enc.splits[2][1]);
#endif
  for (int c=0; c<3; c++) {
    fprintf(stderr,"  bytes %c: %.1f (",
        planesYIQ[c],
        (enc.ctx->bitsSplit[c][0]+enc.ctx->bitsSplit[c][1]+enc.ctx->bitsInterpol[c]+enc.ctx->bitsPixeldata[c][0]+enc.ctx->bitsPixeldata[c][1])*log4k_base/8);
    fprintf(stderr,"%.1f pixels (%.3fb/p)",
        enc.ctx->bitsPixeldata[c][0]*log4k_base/8,
        enc.ctx->symbPixeldata[c][0]>0 ? enc.ctx->bitsPixeldata[c][0]/enc.ctx->symbPixeldata[c][0]*log4k_base : 0.0);
    if (enc.ctx->bitsSplit[c][0]>0)
    fprintf(stderr,", %.1f AS (%.2fb/s)",
        enc.ctx->bitsSplit[c][0]*log4k_base/8,
        enc.ctx->symbSplit[c][0]>0 ? enc.ctx->bitsSplit[c][0]/enc.ctx->symbSplit[c][0]*log4k_base : 0.0);
    if (enc.ctx->bitsSplit[c][1]>0)
    fprintf(stderr,", %.1f LS (%.2fb/s)",
        enc.ctx->bitsSplit[c][1]*log4k_base/8,
        enc.ctx->symbSplit[c][1]>0 ? enc.ctx->bitsSplit[c][1]/enc.ctx->symbSplit[c][1]*log4k_base : 0.0);
    if (enc.ctx->bitsInterpol[c]>0)
    fprintf(stderr,", %.1f IM (%.2fb/s)",
        enc.ctx->bitsInterpol[c]*log4k_base/8,
        enc.ctx->symbInterpol[c]>0 ? enc.ctx->bitsInterpol[c]/enc.ctx->symbInterpol[c]*log4k_base : 0.0);
    if (enc.ctx->bitsPixeldata[c][1]>0)
    fprintf(stderr," , %.1f 1pix (%.2fb/p)",
        enc.ctx->bitsPixeldata[c][1]*log4k_base/8,
        enc.ctx->symbPixeldata[c][1]>0 ? enc.ctx->bitsPixeldata[c][1]/enc.ctx->symbPixeldata[c][1]*log4k_base : 0.0);
    fprintf(stderr,")\n");
  }
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
  fprintf(stderr,"Pixels at depth: \n");
     for(int d=0; d<DEPTH_LEVELS; d++) {
                fprintf(stderr,"depth %i:  %li pixels.\n",d,pixels_at_depth[d]);
     }
*/
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


  if (phase == 1) {
    FILE *g = fopen("context_trees.dat","wb");
    symb_init_write(&enc.coder,g,cutoff);
    write_context_trees(&enc);
    symb_flush(&enc.coder);
    fclose(g);
  }
}

// decode an image
int static decode(image_t *img, char *in) {
  char c[4];
  FILE *f = fopen(in,"rb");
  int ret = fread(c,4,1,f);
  if (ret < 1 || memcmp(c,magic,4) != 0) {
    fprintf(stderr,"No JiF file\n");
    return 3;
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

  int min[3] = {0,0,0};
  int max[3] = {255,510,510};

  ctx_t *ctx=malloc(sizeof(ctx_t));
  coder_t dec = {.source = NULL,.rec = img, .qz={qzy,qzi,qzq}, .method={}, .method_gain={}, .ctx=ctx, .nb_pixels = w*h,
        .range_min={min[0],min[1],min[2]}, .range_max={max[0],max[1],max[2]} };


  dec.minsplitA=0;
  dec.minsplitA |= (fgetc(f) << 0);
  dec.minsplitA |= (fgetc(f) << 8);
  dec.minsplitA |= (fgetc(f) << 16);
  dec.minsplitA |= (fgetc(f) << 24);

  dec.minsplitL=0;
  dec.minsplitL |= (fgetc(f) << 0);
  dec.minsplitL |= (fgetc(f) << 8);
  

  fprintf(stderr,"Decoding %ix%i image with quantization (%i,%i,%i)\n",w,h,qzy/COLOR_STRETCH,qzi/COLOR_STRETCH,qzq/COLOR_STRETCH);
  image_init(img,w,h,&purple);

  dec.sb[0]  |= (fgetc(f) << 0);
  dec.sb[1]  |= (fgetc(f) << 0);
  dec.sb[2]  |= (fgetc(f) << 0);

  if (dec.sb[0] < 8 || dec.sb[1] < 9 || dec.sb[2] < 9) 
          fprintf(stderr,"Significant mantissa bits: Y: %i/8  I: %i/9  Q: %i/9\n", dec.sb[0], dec.sb[1], dec.sb[2]);

  
  if (cutoff != 8)
  fprintf(stderr,"Using RAC range [%g-%g] (cutoff %i)\n",cutoff/4096.0,(4096-cutoff)/4096.0,cutoff);
  symb_init_read(&dec.coder,f,cutoff);
  int phase=1;
  ctx_init(dec.ctx,cutoff,phase);
//  extra_ctx_init(&dec,cutoff);
  input_interpolation_methods(&dec);

  if (dec.minsplitA>0) fprintf(stderr,"Areas smaller than %i pixels are always splitted\n",dec.minsplitA);
  if (dec.minsplitL>0) fprintf(stderr,"Lines shorter than %i pixels are always splitted\n",dec.minsplitL);
  image_set_alpha(img,0);
  int fail = 0;
  for (int i = 0; i < 16; i++) {
    int x = symb_get_simple_bit(&dec.coder);
    fail |= (x != ((imagic>>i) & 1));
  }
  if (fail) fprintf(stderr,"Header corrupt. Will not stop, but result will be bogus.\n");

  fprintf(stderr,"Color ranges: ");
  for(int c=0; c<3; c++) {
        input_min_max(&dec,&(dec.range_min[c]),&(dec.range_max[c]),0,(c==0?255:510),c);
        fprintf(stderr,"%c in [%i,%i]   ",planesYIQ[c],dec.range_min[c],dec.range_max[c]);
  }
  fprintf(stderr,"\n");

  int use_color_data = symb_get_simple_bit(&dec.coder);
  if (use_color_data) {
        dec.color_data.used = 1;
        dec.color_data.max_diff[0] = dec.qz[0]/COLOR_STRETCH;
        dec.color_data.max_diff[1] = dec.qz[1]/COLOR_STRETCH;
        dec.color_data.max_diff[2] = dec.qz[2]/COLOR_STRETCH;
  }


  if (use_color_data) {
        fprintf(stderr,"Reading color table... ");
        input_colors(&dec);
        fprintf(stderr,"done\n");
  }

  int y1 = 0, y2 = img->h-1, x1 = 0, x2 = img->w-1;

  int contexts_given = symb_get_simple_bit(&dec.coder);
  if (contexts_given) {
    fprintf(stderr,"Reading context trees: ");
    read_context_trees(&dec);
    dec.phase = 2;
  }


  fprintf(stderr,"\nHeader read, starting actual decoding.\n");
  progress_bar_firstline();

  int t_m = symb_get_simple_bit(&dec.coder);
  dec.traversal_method = t_m;

  if (t_m == 1) {
    decode_ffv1(&dec,x1,x2,y1,y2);
  } else {
    uint32_t p= img->w * img->h;
    uint32_t bvar[3]={};
    input_pixel(&dec,x1,y1,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
    input_pixel(&dec,x2,y1,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
    input_pixel(&dec,x2,y2,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
    input_pixel(&dec,x1,y2,NULL,p,0,NULL,NULL,NULL,NULL,bvar,bvar,0,0,0);
#if FULL_BORDERS
    int mask=0;
    decode_recurse_line(&dec,x1,x2,y1,y1,mask,0,bvar,bvar,0);
    decode_recurse_line(&dec,x1,x2,y2,y2,mask,0,bvar,bvar,0);
    decode_recurse_line(&dec,x1,x1,y1,y2,mask,0,bvar,bvar,0);
    decode_recurse_line(&dec,x2,x2,y1,y2,mask,0,bvar,bvar,0);
    decode_recurse(&dec,x1+1,x2-1,y1+1,y2-1,0,0,0,0,0,0,SHAPE_0,0,1);
#else
    decode_recurse(&dec,x1,x2,y1,y2,0,0,0,0,0,0,SHAPE_0,0,1);
#endif
  }

  fprintf(stderr,"]\n\n");


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
  yiq2rgb(img);
  fclose(f);
  if (fail & 1) {
    fprintf(stderr,"Image file is invalid\n");
    return 2;
  } else {
    if (fail & 2) {
      fprintf(stderr,"Checksum mismatch\n");
      return 1;
    } else {
      fprintf(stderr,"Checksum matches: 0x%08x\n",crc);
      return 0;
    }
  }
}

void show_help() {
    fprintf(stderr,"Usage: (encoding)\n");
    fprintf(stderr,"   jpif [options] <input.pnm> <output.jpif> [quality_settings]\n");
   
//    fprintf(stderr,"   jpif [-(n)a] [-e] ");
#if FULL_BORDERS    
//    fprintf(stderr,"[-i skipareasplit skiplinesplit]");
#else
//    fprintf(stderr,"[-i skipareasplit]");
#endif
//    fprintf(stderr,"[-c cutoff] [-ffv1] [-m METHODS] [-p PHASE] <input> <output> [quality_y] [quality_i] [quality_q] [factor] [factor_s] [qz_y] [qz_i] [qz_q] [sb_y] [sb_i] [sb_q]\n");
    fprintf(stderr,"Usage: (decoding)\n   jpif -d <input.jpif> <output.pnm>\n");
}

int file_exists(const char * filename){    
        FILE * file = fopen(filename, "r");
        if (file) {        
                fclose(file);        
                return 1;    
        }    
        return 0;
}

// main is just main, right?
int main(int argc, char **argv) {
  char* methods = "2";
  argc--;
  argv++;
  int mode = 0; // 0=encode 1=decode
  int auto_indexing_chart = 0;  // 0=auto,  1=forced yes,  2=forced not
  int minsplitA=-1;
  int minsplitL=0;
  int cutoff=8;
  int tm=0;     // 0=splitting, 1=ffv1
  int phase = 0;        // 0=one-pass   1=phase 1       2=phase 2




  
  while(argc>0) {
          if (strcmp(argv[0],"-e") == 0) {
            mode = 0;
            argc--; argv++; continue; }
          if (strcmp(argv[0],"-d") == 0) {
            mode = 1;
            argc--; argv++; continue; }
          if (strcmp(argv[0],"-a") == 0) {
            auto_indexing_chart = 1;
            argc--; argv++;
            if (argc == 0) {
                fprintf(stderr,"Writing general auto-indexing chart...\n");
                make_auto_indexing_chart();
                return 0;
            }
            continue; }
          if (strcmp(argv[0],"-na") == 0) {
            auto_indexing_chart = 2;
            argc--; argv++; continue; }
          if (strcmp(argv[0],"--help") == 0) {
            show_help(); return 0; }

          if (strcmp(argv[0],"-i") == 0) {
            if (mode == 1) {fprintf(stderr,"Warning: -i option specified while decoding (its value will be ignored).\n");};
            if (tm == 1) {fprintf(stderr,"Warning: -i option specified while using FFV1 traversal (its value will be ignored).\n");};
#if FULL_BORDERS    
            if (argc < 4) {fprintf(stderr,"Option -i expects two numbers\n"); return 1; }
            minsplitA=(int)strtol(argv[1],NULL,10);
            minsplitL=(int)strtol(argv[2],NULL,10);
            argc -= 3;
            argv += 3;
#else
            if (argc < 3) {fprintf(stderr,"Option -i expects a number\n"); return 1; }
            minsplitA=(int)strtol(argv[1],NULL,10);
            argc -= 2;
            argv += 2;
#endif
            continue;            }
          if (strcmp(argv[0],"-c") == 0) {
            if (mode == 1) {fprintf(stderr,"Warning: -c option specified while decoding (its value will be ignored).\n");};
            if (argc < 3) {fprintf(stderr,"Option -c expects a number\n"); return 1; }
            cutoff=(int)strtol(argv[1],NULL,10);
            argc -= 2;
            argv += 2;
            if (cutoff<1 || cutoff>2000) { fprintf(stderr,"Warning: cutoff value should be between 1 and 2000, using default (8)\n"); cutoff=8; }
            continue;          }
          if (strcmp(argv[0],"-ffv1") == 0) {
            if (mode == 1) {fprintf(stderr,"Warning: -ffv1 option specified while decoding (it will be ignored).\n");};
            tm=1;
            argc--; argv++; continue; }
          if (strcmp(argv[0],"-m") == 0) {
            if (mode == 1) {fprintf(stderr,"Warning: -m option specified while decoding (its value will be ignored).\n");};
            if (tm == 1) {fprintf(stderr,"Warning: -m option specified while using FFV1 traversal (its value will be ignored).\n");};
            if (argc < 3) {fprintf(stderr,"Option -m expects a list of methods\n"); return 1;}
            methods = argv[1];
            argc -= 2;
            argv += 2;
            continue;
          }
          if (strcmp(argv[0],"-p") == 0) {
            if (mode == 1) {fprintf(stderr,"Warning: -p option specified while decoding (its value will be ignored).\n");};
            if (argc < 3) {fprintf(stderr,"Option -p expects a number (1 or 2)\n"); return 1; }
            phase=(int)strtol(argv[1],NULL,10);
            argc -= 2; argv += 2;
            if (phase<1 || phase>2) { fprintf(stderr,"Warning: phase should be 1 or 2, using default (1)\n"); phase=1; }
            continue;
          }

          
          if (file_exists(argv[0])) {
                  char *f = strrchr(argv[0],'/');
                  char *ext = f ? strrchr(f,'.') : strrchr(argv[0],'.');
                  if (mode == 0) {
                          if (ext && ( !strcasecmp(ext,".png") ||  !strcasecmp(ext,".pnm") ||  !strcasecmp(ext,".ppm"))) {
                                // ok
                          } else {
                                fprintf(stderr,"Warning: expected \".png\" or \".pnm\" file name extension for input file, trying anyway...\n");
                          }
                  } else {
                          if (ext && ( !strcasecmp(ext,".jpif") ||  !strcasecmp(ext,".jif"))) {
                                // ok
                          } else {
                                fprintf(stderr,"Warning: expected file name extension \".jpif\" for input file, trying anyway...\n");
                          }
                  }
                  break;
          }
          if (argv[0][0] == '-') {
                fprintf(stderr,"Unrecognized option: %s\n",argv[0]);
          } else {
                fprintf(stderr,"Input file does not exist: %s\n",argv[0]);
          }
          show_help();
          return 1;
  }
  if (argc == 0) {
        fprintf(stderr,"Input file missing.\n");
        show_help();
        return 1;
  }
  if (argc == 1) {
        fprintf(stderr,"Output file missing.\n");
        show_help();
        return 1;
  }

  fprintf(stderr,"       ___ __      ____\n");
  fprintf(stderr,"        /  / \\  '  /        |  JPiF 0.1   [May 2011]\n");
  fprintf(stderr,"       /  /__/ /  /-        |  by Jon Sneyers & Pieter Wuille\n");
  fprintf(stderr,"   \\__/  /    /  /          |  (c) 2010-2011;  GNU GPL.\n");
  fprintf(stderr,"____________________________________________________________________\n");

  
  if (mode == 0) {

    // NTSC:  Y: 4 MHz  I: 1.3 MHz  Q: 0.4 MHz

    double q = argc > 2 ? strtod(argv[2],NULL) : 2;
    double qi = argc > 3 ? strtod(argv[3],NULL) : (q<3 ? q*2.5 : 1.5+q*2);
    double qq = argc > 4 ? strtod(argv[4],NULL) : (q<3 ? q*2.5 : 1.5+q*2);
    double factor = argc > 5 ? strtod(argv[5],NULL) : 1.0*(2*q+qi+qq)/6.0;
    double factor_s = argc > 6 ? strtod(argv[6],NULL) : factor/5.0;
    if (factor_s < 1) factor_s = 1;
    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : clamp(1+q/2,1,10),1,100);  // q/2.0-1.0;
    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : clamp(1+qi/2,1,16),1,200); // qi/2.0-1.0;
    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : clamp(1+qq/2,1,18),1,200); // qq/2.0-1.0;
//    double qzy = clamp(argc > 7 ? strtod(argv[7],NULL) : 1+(q+1)/3,1,100);  // q/2.0-1.0;
//    double qzi = clamp(argc > 8 ? strtod(argv[8],NULL) : 1+(qi+1)/3,1,200); // qi/2.0-1.0;
//    double qzq = clamp(argc > 9 ? strtod(argv[9],NULL) : 1+(qq+1)/3,1,200); // qq/2.0-1.0;
    int sb[3] = {};
    sb[0] = clamp(argc > 10 ? strtol(argv[10],NULL,10) : clamp(9 - qzy/2,1,8),8,8);  // using all bits always, otherwise bugs in symbol out
    sb[1] = clamp(argc > 11? strtol(argv[11],NULL,10) : clamp(10 - qzi/2,1,9),9,9);
    sb[2] = clamp(argc > 12? strtol(argv[12],NULL,10) : clamp(10 - qzq/2,1,9),9,9);

    fprintf(stderr,"Using %lu kB of context data\n",(unsigned long) sizeof(ctx_t)/1024);
    fprintf(stderr,"Reading image file: %s\n",argv[0]);
    image_t image;
    if (image_load(argv[0],&image)) {
        fprintf(stderr,"Could not load input image.\n");
        return 1;
    }
    if (minsplitA<0) minsplitA = (q<2 ? image.w*image.h : (q<6 ? 246-40*q : (q<8 ? 12-q : 0)) );
         // clamp((1+64*factor*factor*(2*q+qi+qq)),0,image.w*image.h);

    encode(&image,q,qi,qq,factor,factor_s,qzy,qzi,qzq,argv[1],methods,minsplitA,minsplitL,cutoff,sb,tm,phase,auto_indexing_chart);
    image_free(&image);
//    fprintf(stderr,"maxdists= Y:%.17g  I:%.17g  Q:%.17g\n",maxdist[0],maxdist[1],maxdist[2]);
  } else {
    image_t image;
    int result = decode(&image,argv[0]);
    if (result < 3) image_save(argv[1],&image);
    image_free(&image);
    return result;
  }
  return 0;
}
