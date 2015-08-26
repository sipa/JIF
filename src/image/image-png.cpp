#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <png.h>

#include "image.h"
#include "image-png.h"

int image_load_png(const char *filename, Image &image) {
  FILE *fp = fopen(filename,"rb");
  if (!fp) {
    return 1;
  }
  png_byte header[8];
  int rr = fread(header,1,8,fp);
  int is_png = !png_sig_cmp(header,0,rr);
  if (!is_png) {
    fclose(fp);
    return 2;
  }
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,NULL,NULL);
  if (!png_ptr) {
    fclose(fp);
    return 3;
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_read_struct(&png_ptr,(png_infopp) NULL,(png_infopp) NULL);
    fclose(fp);
    return 4;
  }

  png_init_io(png_ptr,fp);
  png_set_sig_bytes(png_ptr,8);

  png_read_png(png_ptr,info_ptr,PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING | PNG_TRANSFORM_EXPAND, NULL);
//      | PNG_TRANSFORM_STRIP_ALPHA

  size_t width = png_get_image_width(png_ptr,info_ptr);
  size_t height = png_get_image_height(png_ptr,info_ptr);
  int bit_depth = png_get_bit_depth(png_ptr,info_ptr);
  int color_type = png_get_color_type(png_ptr,info_ptr);

  unsigned int nbplanes = 3;
  if (color_type == PNG_COLOR_TYPE_RGB) nbplanes=3;
  else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) nbplanes=4;
  else { printf("Unsupported PNG color type\n"); return 5; }
  image.init(width, height, 0, (1<<bit_depth)-1, nbplanes);

  png_bytepp rows = png_get_rows(png_ptr,info_ptr);

  for (size_t r = 0; r < height; r++) {
    png_bytep row = rows[r];
    for (size_t c = 0; c < width; c++) {
        for (unsigned int p=0; p<nbplanes; p++) {
          image(p,r,c) = (uint8_t) row[c * nbplanes + p];
        }
    }
  }

  png_destroy_read_struct(&png_ptr,&info_ptr,(png_infopp) NULL);
  fclose(fp);

  return 0;
}


int image_save_png(const char *filename, const Image &image) {
  FILE *fp = fopen(filename,"wb");
  if (!fp) {
    return (1);
  }
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,NULL,NULL);
  if (!png_ptr) {
    fclose(fp);
    return (2);
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr,(png_infopp) NULL);
    fclose(fp);
    return (3);
  }

  png_init_io(png_ptr,fp);

  png_set_filter(png_ptr,0,PNG_FILTER_PAETH);
  png_set_compression_level(png_ptr,Z_BEST_COMPRESSION);

  int colortype=PNG_COLOR_TYPE_RGB;
  int nbplanes = image.numPlanes();
  if (nbplanes == 4) colortype=PNG_COLOR_TYPE_RGB_ALPHA;
  if (nbplanes < 3) {printf("FIXME: png with less than 3 planes"); return 4;}
  png_set_IHDR(png_ptr,info_ptr,image.cols(),image.rows(),8,colortype,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,
      PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr,info_ptr);

  png_bytep row = (png_bytep) png_malloc(png_ptr,nbplanes * image.cols());

  for (size_t r = 0; r < (size_t) image.rows(); r++) {
    for (size_t c = 0; c < (size_t) image.cols(); c++) {
      for (int p=0; p<nbplanes; p++) {
        row[c * nbplanes + p] = (png_byte) (image(p,r,c));
      }
    }
    png_write_row(png_ptr,row);
  }

  png_free(png_ptr,row);

  png_write_end(png_ptr,info_ptr);
  png_destroy_write_struct(&png_ptr,&info_ptr);
  fclose(fp);
  return 0;
}
