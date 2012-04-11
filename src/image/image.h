#ifndef _IMAGE_H_
#define _IMAGE_H_ 1

#include <vector>

typedef int ColorVal;

class Plane
{
public:
    unsigned int width, height;
    ColorVal min, max;
    std::vector<ColorVal> data;

    Plane(int width, int height, ColorVal min, ColorVal max) { init(width, height, min, max); }
    Plane() { init(0,0,0,0); }
    void init(int width, int height, ColorVal min, ColorVal max);
    ColorVal &operator()(int r, int c) { return data[r*width + c]; }
    const ColorVal &operator()(int r, int c) const { return data[r*width + c]; }
};

class Image
{
protected:
    std::vector<Plane> planes;
public:

    Image(int width, int height, ColorVal min, ColorVal max, int planes) { init(width, height, min, max, planes); }
    Image() { init(0,0,0,0,0); }

    void init(int width, int height, ColorVal min, ColorVal max, int planes);
    bool load(const char *name);
    bool save(const char *name) const;
    Plane &operator()(int plane) { 
        return planes[plane];
    }
    const Plane &operator()(int plane) const { 
        return planes[plane];
    }
    ColorVal operator()(int plane, int r, int c) const { 
        return planes[plane](r,c);
    }
    ColorVal& operator()(int plane, int r, int c) { 
        return planes[plane](r,c);
    }
    int numPlanes() const { return planes.size(); }
};

#endif