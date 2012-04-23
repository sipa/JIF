#ifndef _IMAGE_H_
#define _IMAGE_H_ 1

#include <vector>
#include <assert.h>

typedef int ColorVal;

class Plane
{
public:
    int subwidth, subheight;
    ColorVal min, max;
    std::vector<ColorVal> data;

    Plane(int subwidth, int subheight, ColorVal min, ColorVal max) {
        init(subwidth, subheight, min, max);
    }

    Plane() {
        init(0,0,0,0);
    }

    void init(int subwidth, int subheight, ColorVal min, ColorVal max);

    ColorVal &operator()(int subr, int subc) {
        assert(!(subr >= subheight || subr < 0 || subc >= subwidth || subc < 0));
        return data[subr*subwidth + subc];
    }

    ColorVal operator()(int subr, int subc) const {
        if (subr >= subheight || subr < 0 || subc >= subwidth || subc < 0) return 0;
        return data[subr*subwidth + subc];
    }
};

class Image
{
protected:
    int width, height;
    std::vector<Plane> planes;
    std::vector<std::pair<int,int> > subsample; // first=rows, seconds=cols

    Plane &operator()(int plane) {
        return planes[plane];
    }
    const Plane &operator()(int plane) const {
        return planes[plane];
    }

public:

    Image(int width, int height, ColorVal min, ColorVal max, int planes) {
        init(width, height, min, max, planes);
    }

    Image() {
        reset();
    }

    void init(int width, int height, ColorVal min, ColorVal max, int planes);

    void reset() {
        init(0,0,0,0,0);
    }

    void add_plane(ColorVal min, ColorVal max, int subSampleR = 1, int subSampleC = 1);

    bool load(const char *name);
    bool save(const char *name) const;

    bool is_set(int p, int r, int c) const {
        return ((r % subsample[p].first) == 0 && (c % subsample[p].second) == 0);
    }

    ColorVal operator()(int p, int r, int c) const {
        return planes[p](r / subsample[p].first,c / subsample[p].second);
    }
    ColorVal& operator()(int p, int r, int c) {
        return planes[p](r / subsample[p].first,c / subsample[p].second);
    }

    int numPlanes() const {
        return planes.size();
    }

    int min(int p) const {
        return planes[p].min;
    }

    int max(int p) const {
        return planes[p].max;
    }

    int rows() const {
        return height;
    }

    int cols() const {
        return width;
    }

    int subSampleR(int p) const {
        return subsample[p].first;
    }

    int subSampleC(int p) const {
        return subsample[p].second;
    }
};

#endif
