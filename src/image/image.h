#ifndef _IMAGE_H_
#define _IMAGE_H_ 1

#include "../config.h"
#include "pixel.h"

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

typedef struct {
  pixel_t *pixels;
  size_t w, h;
} image_t;

void image_free(image_t *image);
void image_set_alpha(image_t *image, uint16_t alpha);
void image_init(image_t *image, int w, int h, const pixel_t *color);

pixel_t static inline *image_pixel(image_t *image, int x, int y) {
  if (x < 0) x = 0;
  if (y < 0) y = 0;
  if (x >= image->w) x = image->w -1;
  if (y >= image->h) y = image->h -1;

  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  
  return image->pixels + (image->w * y + x);
}

xwpixel_t static inline image_pixel_contrast(image_t *image, int x, int y) {
  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  pixel_t *block[9] = {
        image_pixel(image,x-1,y-1), image_pixel(image,x,y-1), image_pixel(image,x+1,y-1), 
        image_pixel(image,x-1,y),   image_pixel(image,x,y),   image_pixel(image,x+1,y), 
        image_pixel(image,x-1,y+1), image_pixel(image,x,y+1), image_pixel(image,x+1,y+1)  
        };
  wpixel_t avg = {};
  for (int i=0; i<9; i++) {
        wpixel_addu(&avg, block[i], 1, 0);
  }
  wpixel_div(&avg, 9, 0);
  xwpixel_t variance = {};
  for (int i=0; i<9; i++) {
        xwpixel_t tmp = {};
        xwpixel_addu(&tmp, block[i], 1, 0);
        xw1wpixel_add(&tmp, &avg, -1, 0);
        xwpixel_square(&tmp, 0);
        xwwpixel_add(&variance, &tmp, 1, 0);
  }
  xwpixel_div(&variance, 9, 0);
  xwpixel_sqrt(&variance, 0);
//  fprintf(stderr,"contrast (%i,%i): %lli, %lli, %lli\n", x,y, variance.wd[0], variance.wd[1], variance.wd[2]);
  return variance; 
}
xwpixel_t static inline image_pixel_gx(image_t *image, int x, int y) {
  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  pixel_t *block[9] = {
        image_pixel(image,x-1,y-1), image_pixel(image,x,y-1), image_pixel(image,x+1,y-1), 
        image_pixel(image,x-1,y),   image_pixel(image,x,y),   image_pixel(image,x+1,y), 
        image_pixel(image,x-1,y+1), image_pixel(image,x,y+1), image_pixel(image,x+1,y+1)  
        };
  xwpixel_t g = {};
  xwpixel_addu(&g, block[0], -1, 0);
  xwpixel_addu(&g, block[1], -2, 0);
  xwpixel_addu(&g, block[2], -1, 0);
  xwpixel_addu(&g, block[3],  2, 0);
  xwpixel_addu(&g, block[4],  4, 0);
  xwpixel_addu(&g, block[5],  2, 0);
  xwpixel_addu(&g, block[6], -1, 0);
  xwpixel_addu(&g, block[7], -2, 0);
  xwpixel_addu(&g, block[8], -1, 0);
  xwpixel_div(&g, 4, 0);
  return g;
}
xwpixel_t static inline image_pixel_gy(image_t *image, int x, int y) {
  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  pixel_t *block[9] = {
        image_pixel(image,x-1,y-1), image_pixel(image,x,y-1), image_pixel(image,x+1,y-1), 
        image_pixel(image,x-1,y),   image_pixel(image,x,y),   image_pixel(image,x+1,y), 
        image_pixel(image,x-1,y+1), image_pixel(image,x,y+1), image_pixel(image,x+1,y+1)  
        };
  xwpixel_t g = {};
  xwpixel_addu(&g, block[0], -1, 0);
  xwpixel_addu(&g, block[1],  2, 0);
  xwpixel_addu(&g, block[2], -1, 0);
  xwpixel_addu(&g, block[3], -2, 0);
  xwpixel_addu(&g, block[4],  4, 0);
  xwpixel_addu(&g, block[5], -2, 0);
  xwpixel_addu(&g, block[6], -1, 0);
  xwpixel_addu(&g, block[7],  2, 0);
  xwpixel_addu(&g, block[8], -1, 0);
  xwpixel_div(&g, 4, 0);
  return g;
}
xwpixel_t static inline image_pixel_nb_distinct(image_t *image, int x, int y) {
  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  pixel_t *block[9] = {
        image_pixel(image,x-1,y-1), image_pixel(image,x,y-1), image_pixel(image,x+1,y-1), 
        image_pixel(image,x-1,y),   image_pixel(image,x,y),   image_pixel(image,x+1,y), 
        image_pixel(image,x-1,y+1), image_pixel(image,x,y+1), image_pixel(image,x+1,y+1)  
        };
  xwpixel_t g = {};
  for (int i=0; i<9; i++) {
    for (int c=0; c<3; c++) {
      int unique = 1;
      for (int j=0; j<i; j++) {
        if (UNSTRETCH(block[i]->d[c]) == UNSTRETCH(block[j]->d[c])) {
                unique = 0;
                break;
        }
      }
      g.wd[c] += unique;
      assert(g.wd[c] <= 10);
    }
  }  
  return g;
}

const pixel_t static inline *image_pixel_c(const image_t *image, int x, int y) {
  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  return image->pixels + (image->w * y + x);
}

int image_load(const char *filename, image_t *image);
int image_save(const char *filename, const image_t *image);

#endif
