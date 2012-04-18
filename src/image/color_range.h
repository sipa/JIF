#ifndef _COLOR_RANGE_H_
#define _COLOR_RANGE_H_ 1

#include <vector>

#include "image.h"

class ColorRanges
{
public:
    virtual ~ColorRanges() {};
    virtual int numPlanes() const =0;
    virtual ColorVal min(int p) const =0;
    virtual ColorVal max(int p) const =0;
    virtual ColorVal min(int p, int r, int c) const { return min(p); }
    virtual ColorVal max(int p, int r, int c) const { return max(p); }
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
    ColorVal min(int p, int r, int c) const { return ranges->min(p,r,c); }
    ColorVal max(int p, int r, int c) const { return ranges->max(p,r,c); }
    bool isStatic() const { return ranges->isStatic(); }
};

const ColorRanges *dupRanges(const ColorRanges *ranges);

#endif
