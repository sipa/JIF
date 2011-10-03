#ifndef _PIXEL_H_
#define _PIXEL_H_ 1

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {
  uint16_t d[3];
  uint16_t a;
} pixel_t;

typedef struct {
  int32_t wd[3];
} wpixel_t;

// force an image's pixel to a particular color (ignores planes in mask), and marking it as done
void static inline pixel_set(pixel_t *out, const pixel_t *pixel, int mask) {
  for (int c = 0; c < 3; c++) {
    if (!(mask & (1 << c))) {
      out->d[c] = pixel->d[c];
    }
  }
}

// force an image's pixel to a particular color (ignores planes in mask), and marking it as done
void static inline wpixel_set(pixel_t *out, const wpixel_t *pixel, int mask) {
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
    if (!(mask & (1 << c))) count->wd[c] += mult*in->d[c];
  }
}

void static inline wwpixel_add(wpixel_t *count, const wpixel_t *in, int mult, int mask) {
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
