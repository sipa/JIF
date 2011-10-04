#ifndef _PIXEL_H_
#define _PIXEL_H_ 1

#include "../config.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "../util.h"

typedef struct {
  uint16_t d[3];
  uint16_t a;
} pixel_t;

typedef struct {
  int32_t wd[3];
} wpixel_t;

typedef struct {
  int64_t wd[3];
} xwpixel_t;

// force an image's pixel to a particular color (ignores planes in mask), and marking it as done
void static inline pixel_set(pixel_t *out, const pixel_t *pixel, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = pixel->d[c];
    }
  }
}

void static inline pixel_median7(pixel_t *out, const pixel_t *a1, const pixel_t *a2, const pixel_t *a3, const pixel_t *a4, const pixel_t *a5, const pixel_t *a6, const pixel_t *a7, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = median7(a1->d[c],a2->d[c],a3->d[c],a4->d[c],a5->d[c],a6->d[c],a7->d[c]);
    }
  }
}
void static inline pixel_median5(pixel_t *out, const pixel_t *a1, const pixel_t *a2, const pixel_t *a3, const pixel_t *a4, const pixel_t *a5, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = median5(a1->d[c],a2->d[c],a3->d[c],a4->d[c],a5->d[c]);
    }
  }
}
void static inline pixel_median3(pixel_t *out, const pixel_t *a1, const pixel_t *a2, const pixel_t *a3, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = median_3(a1->d[c],a2->d[c],a3->d[c]);
    }
  }
}


// force an image's pixel to a particular color (ignores planes in mask), and marking it as done
void static inline wpixel_set(pixel_t *out, const wpixel_t *pixel, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
//      out->d[c] = (pixel->wd[c]<0 ? 0 : pixel->wd[c]);
      out->d[c] = clamp(pixel->wd[c],0,STRETCH((c==0 ? 255 : 510)));
    }
  }
}

void static inline xwpixel_set(pixel_t *out, const xwpixel_t *pixel, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = (pixel->wd[c]<0 ? 0 : pixel->wd[c]);
    }
  }
}


void static inline pixel_linear_w(int f, pixel_t *out, const pixel_t *in1, int d1, const pixel_t *in2, int d2, int mask) {
  assert(f <= 16);
  assert(f >= 0);
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) 
       out->d[c] = (2 *  ((int) (in1->d[c]) * d2 * f + (int) (in2->d[c]) * d1 * (16-f)) + (d1 + d2)*16 ) / (2 * (d1 + d2) * 16);
  }
}


int static inline pixel_equal(const pixel_t *a, const pixel_t *b, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c)) &&  a->d[c] != b->d[c]) return 0;
  }
  return 1;
}

void static inline wpixel_add(wpixel_t *count, const pixel_t *in, int mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        if (count->wd[c] > (1<<30)) fprintf(stderr,"Warning: overflow danger (wpixel_add)\n");
        count->wd[c] += mult*in->d[c];
    }
  }
}
void static inline wpixel_addu(wpixel_t *count, const pixel_t *in, int mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        if (count->wd[c] > (1<<30)) fprintf(stderr,"Warning: overflow danger (wpixel_add)\n");
        count->wd[c] += (c==0?2:1)*mult*UNSTRETCH(in->d[c]);
    }
  }
}
void static inline wpixel_abs(wpixel_t *count, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        count->wd[c] = abs(count->wd[c]);
    }
  }
}

void static inline xwpixel_square(xwpixel_t *count, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        count->wd[c] *= count->wd[c];
        if ((count->wd[c] >> 50) > 0) fprintf(stderr,"Warning: overflow danger (wpixel_square)\n");
    }
  }
}
// set d[c] = max(x,d[c])
void static inline xwpixel_max(xwpixel_t *count, int64_t x, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        if (count->wd[c] < x) count->wd[c] = x;
    }
  }
}
void static inline xwpixel_sqrt(xwpixel_t *count, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        count->wd[c] = sqrt(count->wd[c]);
    }
  }
}
void static inline xwpixel_abs(xwpixel_t *count, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        count->wd[c] = abs(count->wd[c]);
    }
  }
}

void static inline xwpixel_add(xwpixel_t *count, const pixel_t *in, int64_t mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        if ((count->wd[c] >> 50) > 0) fprintf(stderr,"Warning: overflow danger! (xwpixel_add)\n");
        count->wd[c] += mult*in->d[c];
    }
  }
}
void static inline xwpixel_addu(xwpixel_t *count, const pixel_t *in, int64_t mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        if ((count->wd[c] >> 50) > 0) fprintf(stderr,"Warning: overflow danger! (xwpixel_add)\n");
        count->wd[c] += (c==0?2:1)*mult*UNSTRETCH(in->d[c]);
    }
  }
}

void static inline wwpixel_add(wpixel_t *count, const wpixel_t *in, int mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) count->wd[c] += mult*in->wd[c];
  }
}

void static inline xw1wpixel_add(xwpixel_t *count, const wpixel_t *in, int64_t mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) count->wd[c] += mult * ((int64_t) in->wd[c]);
  }
}


void static inline xwwpixel_add(xwpixel_t *count, const xwpixel_t *in, int64_t mult, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) count->wd[c] += mult*in->wd[c];
  }
}

void static inline pixel_div(pixel_t *out, int div, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) out->d[c] = (2 * out->d[c] + div) / (2 * div);
  }
}

void static inline wpixel_div(wpixel_t *out, int div, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) out->wd[c] = (2 * out->wd[c] + div) / (2 * div);
  }
}
void static inline xwwpixel_div(xwpixel_t *out, xwpixel_t *div, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) out->wd[c] = (2 * out->wd[c] + div->wd[c]) / (2 * div->wd[c]);
  }
}
void static inline xwpixel_div(xwpixel_t *out, int64_t div, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) out->wd[c] = (2 * out->wd[c] + div) / (2 * div);
  }
}
void static inline xwpixel_divs(xwpixel_t *out, int64_t div, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) out->wd[c] = (2 * STRETCH(out->wd[c]) + div) / (2 * div);
  }
}
void static inline xwpixel_mul(xwpixel_t *out, int64_t mul, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) out->wd[c] *= mul;
  }
}
void static inline wpixel_clamp(wpixel_t *out, const pixel_t *p1, const pixel_t *p2, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
        int min=(p1->d[c]<p2->d[c]?p1->d[c]:p2->d[c]);
        int max=(p1->d[c]<p2->d[c]?p2->d[c]:p1->d[c]);
        if (out->wd[c] < min) out->wd[c] = min;
        if (out->wd[c] > max) out->wd[c] = max;
    }
  }
}
void pixel_avg(pixel_t *out, const pixel_t **in, const int *weight, int count, int mask);
void pixel_med(pixel_t *out, const pixel_t **in, const int *weight, int count, int mask);

void static inline pixel_quadratic(pixel_t *out,const pixel_t *p0, const pixel_t *p1, const pixel_t *p2, int mask) {
  wpixel_t p={};
  //standard quadratic
/*  wpixel_add(&p,p0,-1,mask);
  wpixel_add(&p,p1,3,mask);
  wpixel_add(&p,p2,1,mask);
  wpixel_div(&p,3,mask);
*/
  wpixel_add(&p,p0,-1,mask);
  wpixel_add(&p,p1,3,mask);
  wpixel_add(&p,p2,1,mask);
  wpixel_div(&p,3,mask);
  
  wpixel_clamp(&p,p1,p2,mask);       // restrict to range between p1 and p2
  wpixel_set(out,&p,mask);
}

void static inline pixel_cubic(pixel_t *out,const pixel_t *p0, const pixel_t *p1, const pixel_t *p2, const pixel_t *p3, int mask) {
  wpixel_t p={};
  wpixel_add(&p,p0,1,mask);
  wpixel_add(&p,p1,-4,mask);
  wpixel_add(&p,p2,6,mask);
  wpixel_add(&p,p3,1,mask);
  wpixel_div(&p,4,mask);
  
  wpixel_clamp(&p,p2,p3,mask);       // restrict to range between p2 and p3
  wpixel_set(out,&p,mask);
}

void static inline pixel_linear(pixel_t *out, const pixel_t *in1, int d1, const pixel_t *in2, int d2, int mask) {
  wpixel_t w={};
  if (d1+d2 == 0) {d1=d2=1;}
  wpixel_add(&w,in1,d2,mask);
  wpixel_add(&w,in2,d1,mask);
  wpixel_div(&w,d1+d2,mask);
  wpixel_set(out,&w,mask);
}

#endif
