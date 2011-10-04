#ifndef _IMAGE_PNM_H_
#define _IMAGE_PNM_H_ 1

#include "image.h"

int image_load_pnm(const char *filename, image_t *image);
int image_save_pnm(const char *filename, const image_t *image);

#endif
