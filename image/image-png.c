#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <png.h>

#include "image.h"
#include "image-png.h"
/*
int image_load_png(const char *filename, image_t *image) {
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

//  png_read_png(png_ptr,info_ptr,PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING | PNG_TRANSFORM_EXPAND
//      | PNG_TRANSFORM_STRIP_ALPHA,NULL);


  png_read_info(png_ptr, info_ptr);

  size_t width = png_get_image_width(png_ptr,info_ptr);
  size_t height = png_get_image_height(png_ptr,info_ptr);
  int bit_depth = png_get_bit_depth(png_ptr,info_ptr);
  int color_type = png_get_color_type(png_ptr,info_ptr);

  image->w = width;
  image->h = height;
  size_t isize=width*height*sizeof(pixel_t);
  image->pixels = malloc(isize);
  fprintf(stderr,"Allocating %llu-byte image: %p\n",(unsigned long long)isize,image->pixels);

  if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_ptr);
  if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_gray_1_2_4_to_8(png_ptr);
  if (bit_depth == 16) png_set_strip_16(png_ptr);
  if (color_type & PNG_COLOR_MASK_ALPHA) png_set_strip_alpha(png_ptr);
  if (bit_depth < 8) png_set_packing(png_ptr);
  if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) png_set_gray_to_rgb(png_ptr);
  

//  png_bytepp rows = png_get_rows(png_ptr,info_ptr);


  png_bytep row = png_malloc(png_ptr,3 * image->w);

  for (size_t r = 0; r < height; r++) {
    png_read_row(png_ptr,row,NULL);
    for (size_t c = 0; c < width; c++) {
      switch (color_type << 8 | bit_depth) {
        case (PNG_COLOR_TYPE_RGB << 8 | 0x08): {
          image->pixels[r * width + c].d[0] = (uint8_t) row[c * 3];
          image->pixels[r * width + c].d[1] = (uint8_t) row[c * 3 + 1];
          image->pixels[r * width + c].d[2] = (uint8_t) row[c * 3 + 2];
          image->pixels[r * width + c].a = 0;
          break;
        }
        case (PNG_COLOR_TYPE_RGB_ALPHA << 8 | 0x08): {
          image->pixels[r * width + c].d[0] = (uint8_t) row[c * 4];
          image->pixels[r * width + c].d[1] = (uint8_t) row[c * 4 + 1];
          image->pixels[r * width + c].d[2] = (uint8_t) row[c * 4 + 2];
          image->pixels[r * width + c].a = (uint8_t) row[c * 4 + 3];
          break;
        }
      }
    }
  }

  png_destroy_read_struct(&png_ptr,&info_ptr,(png_infopp) NULL);
  fclose(fp);

  return 0;
}

*/

int image_load_png(const char *filename, image_t *image) {
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

  png_read_png(png_ptr,info_ptr,PNG_TRANSFORM_STRIP_16 | PNG_TRANSFORM_PACKING | PNG_TRANSFORM_EXPAND
      | PNG_TRANSFORM_STRIP_ALPHA,NULL);

  size_t width = png_get_image_width(png_ptr,info_ptr);
  size_t height = png_get_image_height(png_ptr,info_ptr);
  int bit_depth = png_get_bit_depth(png_ptr,info_ptr);
  int color_type = png_get_color_type(png_ptr,info_ptr);

  image->w = width;
  image->h = height;
  image->pixels = malloc(sizeof(pixel_t) * width * height);

  png_bytepp rows = png_get_rows(png_ptr,info_ptr);

  for (size_t r = 0; r < height; r++) {
    png_bytep row = rows[r];
    for (size_t c = 0; c < width; c++) {
      switch (color_type << 8 | bit_depth) {
        case (PNG_COLOR_TYPE_RGB << 8 | 0x08): {
          image->pixels[r * width + c].d[0] = (uint8_t) row[c * 3];
          image->pixels[r * width + c].d[1] = (uint8_t) row[c * 3 + 1];
          image->pixels[r * width + c].d[2] = (uint8_t) row[c * 3 + 2];
          image->pixels[r * width + c].a = 0;
          break;
        }
        case (PNG_COLOR_TYPE_RGB_ALPHA << 8 | 0x08): {
          image->pixels[r * width + c].d[0] = (uint8_t) row[c * 4];
          image->pixels[r * width + c].d[1] = (uint8_t) row[c * 4 + 1];
          image->pixels[r * width + c].d[2] = (uint8_t) row[c * 4 + 2];
          image->pixels[r * width + c].a = (uint8_t) row[c * 4 + 3];
          break;
        }
      }
    }
  }

  png_destroy_read_struct(&png_ptr,&info_ptr,(png_infopp) NULL);
  fclose(fp);

  return 0;
}


int image_save_png(const char *filename, const image_t *image) {
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

  png_set_IHDR(png_ptr,info_ptr,image->w,image->h,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,
      PNG_FILTER_TYPE_DEFAULT);

  png_write_info(png_ptr,info_ptr);

  png_bytep row = png_malloc(png_ptr,3 * image->w);

  for (size_t r = 0; r < image->h; r++) {
    for (size_t c = 0; c < image->w; c++) {
      row[c * 3] = (png_byte) (image->pixels[r * image->w + c].d[0]);
      row[c * 3 + 1] = (png_byte) (image->pixels[r * image->w + c].d[1]);
      row[c * 3 + 2] = (png_byte) (image->pixels[r * image->w + c].d[2]);
    }
    png_write_row(png_ptr,row);
  }

  png_free(png_ptr,row);

  png_write_end(png_ptr,info_ptr);
  png_destroy_write_struct(&png_ptr,&info_ptr);
  fclose(fp);
  return 0;
}
