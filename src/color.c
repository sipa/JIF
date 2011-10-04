#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "config.h"
#include "image/pixel.h"
#include "image/image.h"
#include "color.h"


// convert in-place from RGB to YIQ
void rgb2yiq(image_t *img) {
  int w = img->w;
  int h = img->h;

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      pixel_rgb2yiq(image_pixel(img,x,y));
    }
  }
}

// convert in-place from YIQ to RGB
void yiq2rgb(image_t *img) {
  int w = img->w;
  int h = img->h;

  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      pixel_yiq2rgb(image_pixel(img,x,y));
    }
  }
}
