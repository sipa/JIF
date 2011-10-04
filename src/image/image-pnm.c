#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <png.h>

#include "image.h"
#include "image-pnm.h"

#define PPMREADBUFLEN 256

int image_load_pnm(const char *filename, image_t *image) {
  FILE *fp = fopen(filename,"rb");
        char buf[PPMREADBUFLEN], *t;
        int r;
  if (!fp) {
    return 1;
  }
  unsigned int width=0,height=0;
  unsigned int maxval=0;
        t = fgets(buf, PPMREADBUFLEN, fp);
        if ( (t == NULL) || ( strncmp(buf, "P6\n", 3) != 0 ) ) {
                fprintf(stderr,"PPM file is not of type P6, cannot read other types.\n");
                fclose(fp);
                return 1;
        }
        do
        { /* Px formats can have # comments after first line */
           t = fgets(buf, PPMREADBUFLEN, fp);
           if ( t == NULL ) return 1;
        } while ( strncmp(buf, "#", 1) == 0 || strncmp(buf, "\n", 1) == 0);
        r = sscanf(buf, "%u %u", &width, &height);
        if ( r < 2 ) {
                fclose(fp);
                return 1;
        }

  char bla;
        r = fscanf(fp, "%u%c", &maxval, &bla);
        if ( (r < 2) || maxval<1 || maxval>255 ) {
                fprintf(stderr,"Invalid PPM file.\n");
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
