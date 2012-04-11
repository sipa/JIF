#include <string.h>
#include <stdio.h>

#include "image.h"
#include "image-pnm.h"

void Plane::init(int width, int height, int min, int max) {
    this->width = width;
    this->height = height;
    this->min = min;
    this->max = max;
    this->data = std::vector<ColorVal>(width*height, 0);
}

bool Image::load(const char *filename) {
  const char *f = strrchr(filename,'/');
  const char *ext = f ? strrchr(f,'.') : strrchr(filename,'.');
//  if (ext && !strcasecmp(ext,".png")) {
//    return image_load_png(filename,*this);
//  }
  if (ext && !strcasecmp(ext,".pnm")) {
    return image_load_pnm(filename,*this);
  }
  if (ext && !strcasecmp(ext,".ppm")) {
    return image_load_pnm(filename,*this);
  }
  fprintf(stderr,"Unknown extension for read from: %s\n",ext ? ext : "(none)");
  return false;
}

bool Image::save(const char *filename) const {
  const char *f = strrchr(filename,'/');
  const char *ext = f ? strrchr(f,'.') : strrchr(filename,'.');
//  if (ext && !strcasecmp(ext,".png")) {
//    return image_save_png(filename,*this);
//  }
  if (ext && !strcasecmp(ext,".pnm")) {
    return image_save_pnm(filename,*this);
  }
  if (ext && !strcasecmp(ext,".ppm")) {
    return image_save_pnm(filename,*this);
  }
  fprintf(stderr,"Unknown extension to write to: %s\n",ext ? ext : "(none)");
  return false;
}

void Image::init(int width, int height, int min, int max, int planes) {
    this->planes = std::vector<Plane>(planes, Plane(width, height, min, max));
}

