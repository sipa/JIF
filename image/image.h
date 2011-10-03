#ifndef _IMAGE_H_
#define _IMAGE_H_ 1

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
  assert(x>=0);
  assert(y>=0);
  assert(x<image->w);
  assert(y<image->h);
  return image->pixels + (image->w * y + x);
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
