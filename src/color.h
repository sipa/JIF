#ifndef _COLOR_H_
#define _COLOR_H_ 1

#include "config.h"

#include "image/image.h"



#define COLOR_DIST(x,y) (abs(UNSTRETCH((x))-UNSTRETCH((y))))

void static inline pixel_rgb2yiq(pixel_t *p) {
  int r = p->d[0];
  int g = p->d[1];
  int b = p->d[2];
  int y = ((r + b) / 2 + g) / 2;
  int i = r - b + 255;
  int q = (r + b) / 2 - g + 255;

  // echte YIQ: y = 0.299*r + 0.587*g + 0.114*b
  //            i = 0.595716*r - 0.274453*g -0.321263*b
  //            q = 0.211456*r - 0.522591*g +0.311135*b

  // benadering: y = 0.25*r + 0.5*g + 0.25*b
  //             i = 1 * r          - 1 * b
  //             q = 0.5*r  - 1*g   + 0.5*b

  p->d[0] = STRETCH(y);
  p->d[1] = STRETCH(i);
  p->d[2] = STRETCH(q);
  
}

void static inline pixel_yiq2rgb(pixel_t *p) {
  int y = UNSTRETCH(p->d[0]);
  int i = UNSTRETCH(p->d[1]);
  int q = UNSTRETCH(p->d[2]);
  int r = y + (q + 2) / 2 + (i + 2) / 2 - 256;
  int g = y - (q + 1) / 2 + 128;
  int b = y + (q + 2) / 2 - (i + 1) / 2;
  if (r < 0) r = 0;
  if (g < 0) g = 0;
  if (b < 0) b = 0;
  if (r > 255) r = 255;
  if (g > 255) g = 255;
  if (b > 255) b = 255;

  p->d[0] = r;
  p->d[1] = g;
  p->d[2] = b;
 
}

void static inline get_range_y(int *min, int *max, int *range_min, int *range_max) {
  *min=range_min[0];
  *max=range_max[0];
  *min=STRETCH(*min);
  *max=STRETCH(*max);
  assert(*min <= *max);
}

void static inline get_range_i(int y, int *min, int *max, int *range_min, int *range_max) {
  int miny,maxy;
  get_range_y(&miny,&maxy, range_min, range_max);
  if (y>maxy) y=maxy;
  if (y<miny) y=miny;

  y=UNSTRETCH(y);
  if (y<63) {
    *min=252-4*y;
    *max=258+4*y;
  } else if (y>=192) {
    *min=3+4*(y-192);
    *max=507-4*(y-192);
  } else {
    *min=0;
    *max=510;
  }
  if (range_min[1] >= *min) *min=range_min[1];
  if (range_max[1] <= *max) *max=range_max[1];
  if (range_min[1] >= *max) *max=range_min[1];
  if (range_max[1] <= *min) *min=range_max[1];
  *min=STRETCH(*min);
  *max=STRETCH(*max);
  assert(*min <= *max);
}


void static inline get_range_q(int y, int i, int *min, int *max, int *range_min, int *range_max) {
  int miny,maxy;
  get_range_y(&miny,&maxy, range_min, range_max);
  if (y>maxy) y=maxy;
  if (y<miny) y=miny;

  int mini,maxi;
  get_range_i(y,&mini,&maxi, range_min, range_max);
  if (i>maxi) i=maxi;
  if (i<mini) i=mini;

  y=UNSTRETCH(y);
  i=UNSTRETCH(i);
  if (y<63) {
    *min=254-2*y+(abs(i-255)/2)*2;
    *max=256+2*y;
  } else if (y>=192) {
    *min=255-2*(255-y);
    *max=255+2*(255-y)-((1+abs(i-255))/2)*2;
  } else {
    *min=MAX(1+(y-128)*2,128-(y-63)*2+(abs(i-255)/2)*2);
    *max=MIN(382+(y-63)*2,383+(191-y)*2-((1+abs(i-255))/2)*2);
  }
  if (range_min[2] >= *min) *min=range_min[2];
  if (range_max[2] <= *max) *max=range_max[2];
  if (range_min[2] >= *max) *max=range_min[2];
  if (range_max[2] <= *min) *min=range_max[2];
  *min=STRETCH(*min);
  *max=STRETCH(*max);
  assert(*min <= *max);
}

void rgb2yiq(image_t *img);
void yiq2rgb(image_t *img);

#endif
