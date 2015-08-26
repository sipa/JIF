#ifndef _COLOR_RANGE_H_
#define _COLOR_RANGE_H_ 1

#include <vector>

#include "image.h"

typedef std::vector<ColorVal> prevPlanes;


class ColorRanges
{
public:
    virtual ~ColorRanges() {};
    virtual int numPlanes() const =0;
    virtual ColorVal min(int p) const =0;
    virtual ColorVal max(int p) const =0;
    virtual void minmax(const int p, const prevPlanes &pp, ColorVal &minv, ColorVal &maxv) const { minv=min(p); maxv=max(p); }
    virtual void snap(const int p, const prevPlanes &pp, ColorVal &minv, ColorVal &maxv, ColorVal &v) const {
        minmax(p,pp,minv,maxv);
        if(v>maxv) v=maxv;
        if(v<minv) v=minv;
    }
    virtual ColorVal reconstruct(const int p, const prevPlanes &pp, ColorVal v) const {
        return v;
    }
    virtual ColorVal position(const int p, const prevPlanes &pp, ColorVal v) const {
        return v;
    }

/*
    virtual ColorVal snap(int p, ColorVal v) const { if (v>max(p)) return max(p); if (v<min(p)) return min(p); return v;};
    virtual ColorVal snap(int p, ColorVal v, int r, int c) const { if (v>max(p,r,c)) return max(p,r,c); if (v<min(p,r,c)) return min(p,r,c); return v;};
    virtual ColorVal snap(int p, ColorVal v, int z, int r, int c) const { if (v>max(p,z,r,c)) return max(p,z,r,c); if (v<min(p,z,r,c)) return min(p,z,r,c); return v;};
    virtual ColorVal min(int p, int r, int c) const { return min(p); }
    virtual ColorVal max(int p, int r, int c) const { return max(p); }
    virtual ColorVal min(int p, int z, int r, int c) const { return min(p); }
    virtual ColorVal max(int p, int z, int r, int c) const { return max(p); }
*/
    virtual bool isStatic() const { return true; }
};

typedef std::vector<std::pair<ColorVal, ColorVal> > StaticColorRangeList;

class StaticColorRanges : public ColorRanges
{
protected:
    StaticColorRangeList ranges;

public:
    StaticColorRanges(StaticColorRangeList ranges) { this->ranges = ranges; }
    int numPlanes() const { return ranges.size(); }
    ColorVal min(int p) const { return ranges[p].first; }
    ColorVal max(int p) const { return ranges[p].second; }
};

const ColorRanges *getRanges(const Image &image);

class DupColorRanges : public ColorRanges {
protected:
    const ColorRanges *ranges;
public:
    DupColorRanges(const ColorRanges *rangesIn) : ranges(rangesIn) {}

    int numPlanes() const { return ranges->numPlanes(); }
    ColorVal min(int p) const { return ranges->min(p); }
    ColorVal max(int p) const { return ranges->max(p); }
    void minmax(const int p, const prevPlanes &pp, ColorVal &minv, ColorVal &maxv) const { ranges->minmax(p,pp,minv,maxv); }
//    ColorVal min(int p, int r, int c) const { return ranges->min(p,r,c); }
//    ColorVal max(int p, int r, int c) const { return ranges->max(p,r,c); }
    bool isStatic() const { return ranges->isStatic(); }
};

const ColorRanges *dupRanges(const ColorRanges *ranges);

#endif
