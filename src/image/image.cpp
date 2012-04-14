#include <string.h>
#include <stdio.h>

#include "image.h"
#include "image-pnm.h"

void Plane::init(int subwidth, int subheight, int min, int max)
{
    this->subwidth = subwidth;
    this->subheight = subheight;
    this->min = min;
    this->max = max;
    this->data = std::vector<ColorVal>(subwidth*subheight, 0);
}

bool Image::load(const char *filename)
{
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

bool Image::save(const char *filename) const
{
    const char *f = strrchr(filename,'/');
    const char *ext = f ? strrchr(f,'.') : strrchr(filename,'.');
//  if (ext && !strcasecmp(ext,".png")) {
//    return image_save_png(filename,*this);
//  }
    if (ext && !strcasecmp(ext,".pnm")) {
        return image_save_pnm(filename,*this);
    }
    if (ext && !strcasecmp(ext,".pgm")) {
        return image_save_pnm(filename,*this);
    }
    if (ext && !strcasecmp(ext,".ppm")) {
        return image_save_pnm(filename,*this);
    }
    fprintf(stderr,"Unknown extension to write to: %s\n",ext ? ext : "(none)");
    return false;
}

void Image::add_plane(int min, int max, int subSampleR, int subSampleC)
{
    planes.push_back(Plane((width+subSampleR-1)/subSampleR, (height+subSampleC-1)/subSampleC, min, max));
    subsample.push_back(std::make_pair(subSampleR, subSampleC));
}

void Image::init(int width, int height, int min, int max, int planes)
{
    this->planes = std::vector<Plane>(planes, Plane(width, height, min, max));
    this->width = width;
    this->height = height;
    this->subsample = std::vector<std::pair<int,int> >(planes, std::make_pair(1,1));
}

