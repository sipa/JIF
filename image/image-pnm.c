#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <png.h>

#include "image.h"
#include "image-pnm.h"

int image_load_pnm(const char *filename, image_t *image) {
  FILE *fp = fopen(filename,"rb");
  if (!fp) {
    return 1;
  }
  unsigned int width=0,height=0;
  unsigned int maxval=0;
  char bla;
  int ret=fscanf(fp,"P6 %u %u %u%c",&width,&height,&maxval,&bla);
  if (ret!=4 || maxval<1 || maxval>255) {
    fclose(fp);
    return 2;
  }
  image->w=width;
  image->h=height;
  image->pixels=malloc(sizeof(pixel_t)*width*height);
  for (int y=0; y<image->h; y++) {
    for (int x=0; x<image->w; x++) {
      for (int c=0; c<3; c++) {
        image_pixel(image,x,y)->d[c]=(fgetc(fp)*255+(maxval/2))/maxval;
      }
    }
  }
  fclose(fp);
  return 0;
}

int image_save_pnm(const char *filename, const image_t *image) {
  FILE *fp = fopen(filename,"wb");
  if (!fp) {
    return (1);
  }
  fprintf(fp,"P6\n%lu %lu\n255\n",(unsigned long) (image->w),(unsigned long) (image->h));
  for (int y = 0; y < image->h; y++) {
    for (int x = 0; x < image->w; x++) {
      fputc(image_pixel_c(image,x,y)->d[0] & 0xFF,fp);
      fputc(image_pixel_c(image,x,y)->d[1] & 0xFF,fp);
      fputc(image_pixel_c(image,x,y)->d[2] & 0xFF,fp);
    }
  }
  fclose(fp);
  return 0;
}
