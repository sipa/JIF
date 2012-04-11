#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "image.h"
#include "image-pnm.h"

#define PPMREADBUFLEN 256

bool image_load_pnm(const char *filename, Image& image) {
  FILE *fp = fopen(filename,"rb");
        char buf[PPMREADBUFLEN], *t;
        int r;
  if (!fp) {
    return false;
  }
  unsigned int width=0,height=0;
  unsigned int maxval=0;
        t = fgets(buf, PPMREADBUFLEN, fp);
        if ( (t == NULL) || ( strncmp(buf, "P6\n", 3) != 0 ) ) {
                fprintf(stderr,"PPM file is not of type P6, cannot read other types.\n");
                fclose(fp);
                return false;
        }
        do
        { /* Px formats can have # comments after first line */
           t = fgets(buf, PPMREADBUFLEN, fp);
           if ( t == NULL ) return 1;
        } while ( strncmp(buf, "#", 1) == 0 || strncmp(buf, "\n", 1) == 0);
        r = sscanf(buf, "%u %u", &width, &height);
        if ( r < 2 ) {
                fclose(fp);
                return false;
        }

  char bla;
        r = fscanf(fp, "%u%c", &maxval, &bla);
        if ( (r < 2) || maxval<1 || maxval>255 ) {
                fprintf(stderr,"Invalid PPM file.\n");
                fclose(fp);
                return 2;
        }

  image.init(width, height, 0, maxval, 3);

  for (unsigned int y=0; y<height; y++) {
    for (unsigned int x=0; x<width; x++) {
      for (unsigned int c=0; c<3; c++) {
        image(c,y,x) = fgetc(fp);
      }
    }
  }
  fclose(fp);
  return true;
}

bool image_save_pnm(const char *filename, const Image& image) {
  FILE *fp = fopen(filename,"wb");
  if (!fp) {
    return false;
  }
  
  if (image.numPlanes() != 3)
  {
      fprintf(stderr,"Cannot store as PPM. Find out why.\n");
      fclose(fp);
      return false;
  }

  ColorVal max = std::max(std::max(image(0).max, image(1).max), image(2).max);
  ColorVal min = std::min(std::min(image(0).min, image(1).min), image(2).min);
  
  if (max-min > 255)
  {
      fprintf(stderr,"Cannot store as PPM. Find out why.\n");
      fclose(fp);
      return false;
  }
  
  unsigned int height = image(0).height, width = image(0).width;
  fprintf(fp,"P6\n%u %u\n%i\n", width, height, max-min);
  for (unsigned int y = 0; y < height; y++) {
    for (unsigned int x = 0; x < width; x++) {
      fputc((image(0,y,x) - min) & 0xFF,fp);
      fputc((image(1,y,x) - min) & 0xFF,fp);
      fputc((image(2,y,x) - min) & 0xFF,fp);
    }
  }
  fclose(fp);
  return true;
}
