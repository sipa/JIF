#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <strings.h>

#include "image.h"

#include "image-png.h"
#include "image-pnm.h"

void image_free(image_t *image) {
  free(image->pixels);
  image->w = 0;
  image->h = 0;
  image->pixels = NULL;
}

void image_set_alpha(image_t *image, uint16_t alpha) {
  for (int y = 0; y < image->h; y++) {
    for (int x = 0; x < image->w; x++) {
      image_pixel(image,x,y)->a = alpha;
    }
  }
}

void image_init(image_t *image, int w, int h, const pixel_t *color) {
  image->w = w;
  image->h = h;
  image->pixels = malloc(sizeof(pixel_t) * w * h);
  for (int y = 0; y < image->h; y++) {
    for (int x = 0; x < image->w; x++) {
      *(image_pixel(image,x,y)) = (*color);
    }
  }
}

int image_load(const char *filename, image_t *image) {
  char *f = strrchr(filename,'/');
  char *ext = f ? strrchr(f,'.') : strrchr(filename,'.');
  if (ext && !strcasecmp(ext,".png")) {
    return image_load_png(filename,image);
  }
  if (ext && !strcasecmp(ext,".pnm")) {
    return image_load_pnm(filename,image);
  }
  if (ext && !strcasecmp(ext,".ppm")) {
    return image_load_pnm(filename,image);
  }
  fprintf(stderr,"Unknown extension for read from: %s\n",ext ? ext : "(none)");
  return 3;
}

int image_save(const char *filename, const image_t *image) {
  char *f = strrchr(filename,'/');
  char *ext = f ? strrchr(f,'.') : strrchr(filename,'.');
  if (ext && !strcasecmp(ext,".png")) {
    return image_save_png(filename,image);
  }
  if (ext && !strcasecmp(ext,".pnm")) {
    return image_save_pnm(filename,image);
  }
  if (ext && !strcasecmp(ext,".ppm")) {
    return image_save_pnm(filename,image);
  }
  fprintf(stderr,"Unknown extension to write to: %s\n",ext ? ext : "(none)");
  return 3;
}
