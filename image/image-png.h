#ifndef _IMAGE_PNG_H_
#define _IMAGE_PNG_H_ 1

#include "image.h"

int image_load_png(const char *filename, image_t *image);
int image_save_png(const char *filename, const image_t *image);

#endif
