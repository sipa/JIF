#ifndef _CB_H_
#define _CB_H_ 1

#include <vector>

#include "transform.h"
#include "../maniac/symbol.h"

#define MAX_PER_BUCKET_0 255
#define MAX_PER_BUCKET_1 510
// if bucket0 and bucket1 are big enough to be exact (and no simplification happens!), then we don't need the bit to check for empty buckets
//#define MAX_PER_BUCKET_1 20

#define MAX_PER_BUCKET_2 5
#define MAX_PER_BUCKET_3 255
const unsigned int max_per_colorbucket[] = {MAX_PER_BUCKET_0, MAX_PER_BUCKET_1, MAX_PER_BUCKET_2, MAX_PER_BUCKET_3};

// quantization constants
#define CB0a 1
#define CB0b 1
#define CB1 4


static int totaldiscretecolors=0;
static int totalcontinuousbuckets=0;

class ColorBucket {
public:
    ColorVal min;
    ColorVal max;
    std::vector<ColorVal> values;
    bool discrete;

    ColorBucket() {
        min = 10000;  // +infinity
        max = -1;  // -infinity    (set to empty interval to start with)
        discrete = true;
    }
    void addColor(const ColorVal c, const unsigned int max_per_bucket) {
        if (c<min) min=c;
        if (c>max) max=c;
        if (discrete) {
          unsigned int pos=0;
          for(; pos < values.size(); pos++) {
                if (c == values[pos]) return;
                if (values[pos] > c) break;
          }
          if (values.size() < max_per_bucket) {
                values.insert(values.begin()+pos, c);
                totaldiscretecolors++;
          } else {
                totaldiscretecolors -= max_per_bucket;
                values.clear();
                discrete=false;
                totalcontinuousbuckets++;
          }
        }
    }

    void removeColor(const ColorVal c) {
        if (discrete) {
          unsigned int pos=0;
          for(; pos < values.size(); pos++) {
                if (c == values[pos]) {
                   values.erase(values.begin()+pos);
                   break;
                }
          }
          if (c==min) min = values[0];
          if (c==max) max = values[values.size()-1];
        } else {
          if (c==min) min++;
          if (c==max) max--;
          if (c>max) return;
          if (c<min) return;
          discrete=true;
          values.clear();
          for (ColorVal x=min; x <= max; x++) {
            if (x != c) values.push_back(x);
          }
        }
    }

    bool empty() const {
       return (min>max);
    }
    void simplify_lossless() {
        if (discrete) {
                if ((int)values.size() == max-min+1) {
                        discrete=false;  // bucket actually contains a continuous range
                        totaldiscretecolors -= values.size();
                        totalcontinuousbuckets++;
                        values.clear();
                }
        }
    }
    void simplify(int percent) {
        if (empty()) return;
        simplify_lossless();
        if (discrete) {
                // heuristic: turn discrete bucket into a continuous one if it is dense enough
                if ((int)values.size()-2 > (max-min-1)*percent/100) {
                        discrete=false; // more than <percent> of the ]min,max[ values are present
                        totaldiscretecolors -= values.size();
                        totalcontinuousbuckets++;
                        values.clear();
                }
        }
    }

    ColorVal snapColor(const ColorVal c) const {
        if (c<min) return min;
        if (c>max) return max;
        if (discrete) {
          ColorVal mindiff = abs(c-min);
          unsigned int best = 0;
          for(unsigned int i=0; i < values.size(); i++) {
                if (c == values[i]) return c;
                ColorVal diff = abs(c-values[i]);
                if (diff < mindiff) {best = i; mindiff = diff;}
                if (values[i] > c) break; // can safely skip the rest, values is sorted
          }
          return values[best];
        }
        return c;
    }
/*
    ColorVal lowestIndex() const {
        if (discrete) return 0;
        else return min;
    }
    ColorVal highestIndex() const {
        if (discrete) return values.size()-1;
        else return max;
    }
    ColorVal colorAtIndex(ColorVal c) const {
        if (discrete) return values[c];
        else return c;
    }
    ColorVal findColor(const ColorVal c) const {
        if (discrete) {
          if (c<min) return 0;
          if (c>max) return highestIndex();
          ColorVal mindiff = abs(c-min);
          ColorVal best = 0;
          for(unsigned int i=0; i < values.size(); i++) {
                if (c == values[i]) return i;
                ColorVal diff = abs(c-values[i]);
                if (diff < mindiff) {best = i; mindiff = diff;}
                if (values[i] > c) break; // can safely skip the rest, values is sorted
          }
          return best;
        }
        if (c<min) return min;
        if (c>max) return max;
        return c;
    }
*/
    void print() const {
        if (min>max) printf("E ");
        else if (min==max) printf("S%i ",min);
        else {
                if (discrete) {
                        printf("D[");
                        for (ColorVal c : values) printf("%i%c",c,(c<max?',':']'));
                        printf(" ");
                } else {
                        printf("C%i..%i ",min,max);
                }
        }
    }
    void printshort() const {
        if (min>max) printf(".");
        else if (min==max) printf("1");
        else {
                if (discrete) {
                        if (values.size()>9) printf("+"); else printf("%lu",values.size());
                } else {
                        printf("C");
                }
        }
    }
};

class ColorBuckets {
public:
    ColorBucket bucket0;
    int min0, min1;
    std::vector<ColorBucket> bucket1;
    std::vector<std::vector<ColorBucket> > bucket2;
    ColorBucket bucket3;
    const ColorRanges *ranges;
    ColorBuckets(const ColorRanges *r) : bucket0(), min0(r->min(0)), min1(r->min(1)),
                                         bucket1((r->max(0) - min0)/CB0a +1),
                                         bucket2((r->max(0) - min0)/CB0b +1, std::vector<ColorBucket>((r->max(1) - min1)/CB1 +1)),
                                         bucket3(),
                                         ranges(r) {}
    ColorBucket& findBucket(const int p, const prevPlanes &pp) {
        if (p==0) return bucket0;
        if (p==1) return bucket1[(pp[0]-min0)/CB0a];
        if (p==2) return bucket2[(pp[0]-min0)/CB0b][(pp[1]-min1)/CB1];
        if (p==3) return bucket3;
        printf("OOPS: fixme: %i\n",p);
        return bucket0;
    }
    ColorBucket findBucket(const int p, const prevPlanes &pp) const {
        if (p==0) return bucket0;
        if (p==1) return bucket1[(pp[0]-min0)/CB0a];
        if (p==2) return bucket2[(pp[0]-min0)/CB0b][(pp[1]-min1)/CB1];
        if (p==3) return bucket3;
        printf("OOPS: fixmeB\n");
        return bucket0;
    }
    void addColor(const std::vector<ColorVal> &pixel) {
        for (unsigned int p=0; p < pixel.size(); p++) {
                findBucket(p, pixel).addColor(pixel[p],max_per_colorbucket[p]);
        }
    }

    bool exists(const int p, const prevPlanes &pp) const {
        if (p>0 && (pp[0] < min0 || pp[0] > ranges->max(0))) return false;
        if (p>1 && (pp[1] < min1 || pp[1] > ranges->max(1))) return false;

        ColorVal rmin, rmax;
        ColorVal v=pp[p];
//        printf("check exist pp="); for (auto x : pp) printf("%i ",x); printf("\n");
        ranges->snap(p,pp,rmin,rmax,v);
//        printf("result=%i",v == pp[p]);
        if (v != pp[p]) return false;   // bucket empty because of original range constraints

        const ColorBucket b = findBucket(p,pp);
        //if (b.min > b.max) return false;
        if (b.snapColor(pp[p]) != pp[p]) return false;
        return true;
    }
    bool exists(const int p, const prevPlanes &lower, const prevPlanes &upper) const {
        prevPlanes pixel=lower;
        if (p==0) {
           for (pixel[0]=lower[0]; pixel[0] <= upper[0]; pixel[0]++)
                if (exists(p,pixel)) return true;
        }
        if (p==1) {
           for (pixel[0]=lower[0]; pixel[0] <= upper[0]; pixel[0]++) {
             for (pixel[1]=lower[1]; pixel[1] <= upper[1]; pixel[1]++) {
                if (exists(p,pixel)) return true;
             }
           }
        }
        return false;
    }
    void print() {
        printf("Y buckets:\n");
        bucket0.print();
        printf("\nI buckets:\n");
        for (auto b : bucket1) b.print();
        printf("\nQ buckets:\n  ");
        for (auto bs : bucket2) { for (auto b : bs) b.print(); printf("\n  ");}
        if (ranges->numPlanes() > 3) {
          printf("Alpha buckets:\n");
          bucket3.print();
        }
    }
};

class ColorRangesCB : public ColorRanges
{
protected:
    const ColorRanges *ranges;
    ColorBuckets *buckets;
public:
    ColorBucket bucket(const int p, const prevPlanes &pp) const { return buckets->findBucket(p,pp); }
//    ColorBucket& bucket(const int p, const prevPlanes &pp) { return buckets->findBucket(p,pp); }
    ColorRangesCB(const ColorRanges *rangesIn, ColorBuckets *cbIn) :  ranges(rangesIn), buckets(cbIn) {} //print();}
    bool isStatic() const { return false; }
    int numPlanes() const { return ranges->numPlanes(); }
    ColorVal min(int p) const { return ranges->min(p); }
    ColorVal max(int p) const { return ranges->max(p); }
/*
    ColorVal reconstruct(const int p, const prevPlanes &pp, ColorVal v) const {
        return bucket(p,pp).colorAtIndex(v);
    }
    ColorVal position(const int p, const prevPlanes &pp, ColorVal v) const {
        return bucket(p,pp).findColor(v);
    }
*/
    void snap(const int p, const prevPlanes &pp, ColorVal &minv, ColorVal &maxv, ColorVal &v) const {
        const ColorBucket& b = bucket(p,pp);
/*
        minv=b.lowestIndex();
        maxv=b.highestIndex();
        v=b.findColor(v);
*/
        minv=b.min;
        maxv=b.max;
        if (b.min > b.max) { printf("UGH!! HOW?\n"); assert(false); minv=0; maxv=0; v=0; }
        v=b.snapColor(v);

//        while (true) {
//          v=b.snapColor(v);
/*          if (p==1) {
            prevPlanes npp(pp[0],v);
            if (bucket(2,npp).empty()) b.removeColor(v);
            else break;
          } else break;
        }
*/
    }
    void minmax(const int p, const prevPlanes &pp, ColorVal &minv, ColorVal &maxv) const {
        const ColorBucket& b = bucket(p,pp);
        minv=b.min;
        maxv=b.max;
        assert(minv <= maxv);
    }
    void print() {
        buckets->print();
    }
};


class TransformCB : public Transform {
protected:
    ColorBuckets *cb;

    const ColorRanges* meta(Image& image, const ColorRanges *srcRanges) {
//        cb->print();
        // in the I buckets, some discrete buckets may have become continuous to keep the colorbucket info small
        // this means some Q buckets are empty, which means that some values from the I buckets can be eliminated
        prevPlanes pixelL,pixelU;
        pixelL.push_back(cb->min0);
        pixelU.push_back(cb->min0+CB0b-1);
        pixelL.push_back(cb->min1);
        pixelU.push_back(cb->min1+CB1-1);
        for (auto bv : cb->bucket2) {
                pixelL[1] = cb->min1;
                pixelU[1] = cb->min1+CB1-1;
                for (auto b : bv) {
                        if (b.empty()) {
//                                printf("tweaking bucket\n");
                                for (ColorVal c=pixelL[1]; c<=pixelU[1]; c++) {
                                  cb->findBucket(1,pixelL).removeColor(c);
                                  cb->findBucket(1,pixelU).removeColor(c);
                                }
                        }
                        pixelL[1] += CB1; pixelU[1] += CB1;
                }
                pixelL[0] += CB0b; pixelU[0] += CB0b;
        }

//        cb->print();
        return new ColorRangesCB(srcRanges, cb);
    }
    bool init(const ColorRanges *srcRanges) { if(srcRanges->numPlanes() < 3) return false; cb = new ColorBuckets(srcRanges); return true; }

    void minmax(const ColorRanges *srcRanges, const int p, const prevPlanes &lower, const prevPlanes &upper, ColorVal &smin, ColorVal &smax) const {
        ColorVal rmin, rmax;
        smin=10000;
        smax=-10000;
        prevPlanes pixel=lower;
        if (p==0) {
                srcRanges->minmax(p,pixel,smin,smax);
        } else if (p==1) {
           for (pixel[0]=lower[0]; pixel[0] <= upper[0]; pixel[0]++) {
                srcRanges->minmax(p,pixel,rmin,rmax);
                if (rmin<smin) smin=rmin;
                if (rmax>smax) smax=rmax;
           }
        } else if (p==2) {
           for (pixel[0]=lower[0]; pixel[0] <= upper[0]; pixel[0]++) {
             for (pixel[1]=lower[1]; pixel[1] <= upper[1]; pixel[1]++) {
                srcRanges->minmax(p,pixel,rmin,rmax);
                if (rmin<smin) smin=rmin;
                if (rmax>smax) smax=rmax;
             }
           }
        } else if (p==3) {
                srcRanges->minmax(p,pixel,smin,smax);
        }
    }

    ColorBucket load_bucket(RacIn &rac, const ColorRanges *srcRanges, const int plane, const prevPlanes &pixelL, const prevPlanes &pixelU) {
        ColorBucket b;
        if (plane<3)
        for (int p=0; p<plane; p++) {
                if (!cb->exists(p,pixelL,pixelU)) return b;
        }
        SimpleSymbolCoder<StaticBitChance, RacIn, 24> coder(rac);

        ColorVal smin,smax;
        minmax(srcRanges,plane,pixelL,pixelU,smin,smax);
        if (smin == smax) {b.min = b.max = smin; b.discrete=false; return b;}

        int exists = coder.read_int(0, 1);
        if (exists == 0) {
           return b; // empty bucket
        }
        b.min = coder.read_int(smin, smax);
        b.max = coder.read_int(b.min, smax);
        if (b.min == b.max) { b.discrete=false; return b; }
        if (b.min + 1 == b.max) { b.discrete=false; return b; }
        b.discrete = coder.read_int(0,1);
        if (b.discrete) {
           int nb = coder.read_int(2, max_per_colorbucket[plane]);
           b.values.push_back(b.min);
           ColorVal v=b.min;
           for (int p=1; p < nb-1; p++) {
             b.values.push_back(coder.read_int(v+1, b.max-1));
//             printf("val %i\n",b.values[p]);
             v = b.values[p];
           }
           if (b.min < b.max) b.values.push_back(b.max);
        }
        return b;
    }
    void load(const ColorRanges *srcRanges, RacIn &rac) {
        printf("Loading Color Buckets\n");
        prevPlanes pixelL, pixelU;
        cb->bucket0 = load_bucket(rac, srcRanges, 0, pixelL, pixelU);
        pixelL.push_back(cb->min0);
        pixelU.push_back(cb->min0+CB0a-1);
        for (ColorBucket& b : cb->bucket1) {b=load_bucket(rac, srcRanges, 1, pixelL, pixelU); pixelL[0] += CB0a; pixelU[0] += CB0a; }
        pixelL[0] = cb->min0;
        pixelU[0] = cb->min0+CB0b-1;
        pixelL.push_back(cb->min1);
        pixelU.push_back(cb->min1+CB1-1);
        for (auto& bv : cb->bucket2) {
                pixelL[1] = cb->min1;
                pixelU[1] = cb->min1+CB1-1;
                for (ColorBucket& b : bv) {
                        b=load_bucket(rac, srcRanges, 2, pixelL, pixelU);
                        pixelL[1] += CB1; pixelU[1] += CB1; 
                }
                pixelL[0] += CB0b; pixelU[0] += CB0b; 
        }
        if (srcRanges->numPlanes() > 3) cb->bucket3 = load_bucket(rac, srcRanges, 3, pixelL, pixelU);
    }

    void save_bucket(const ColorBucket &b, RacOut &rac, const ColorRanges *srcRanges, const int plane, const prevPlanes &pixelL, const prevPlanes &pixelU) const {
        if (plane<3)
        for (int p=0; p<plane; p++) {
                if (!cb->exists(p,pixelL,pixelU)) {
                        if (!b.empty()) {printf("\nBucket does not exist but is not empty!\n"); assert(false);}
                        return;
                }
        }
//        b.print();
        SimpleSymbolCoder<StaticBitChance, RacOut, 24> coder(rac);
//        if (b.min > b.max) printf("SHOULD NOT HAPPEN!\n");

        ColorVal smin,smax;
        minmax(srcRanges,plane,pixelL,pixelU,smin,smax);
        if (smin==smax) { return;}

        b.printshort();

        if (b.min > b.max) {
                coder.write_int(0, 1, 0);  // empty bucket
                return;
        } else coder.write_int(0, 1, 1);  // non-empty bucket


        coder.write_int(smin, smax, b.min);
        coder.write_int(b.min, smax, b.max);
        if (b.min == b.max) return;  // singleton bucket
        if (b.min + 1  == b.max) return; // bucket contains two consecutive values
        coder.write_int(0, 1, b.discrete);
        if (b.discrete) {
           coder.write_int(2, max_per_colorbucket[plane], b.values.size());
           ColorVal v=b.min;
           for (unsigned int p=1; p < b.values.size() - 1; p++) {
             coder.write_int(v+1, b.max-1, b.values[p]);
             v = b.values[p];
           }
        }
    }
    void save(const ColorRanges *srcRanges, RacOut &rac) const {
        printf("Saving Y Color Bucket: ");
        prevPlanes pixelL, pixelU;
        save_bucket(cb->bucket0, rac, srcRanges, 0, pixelL, pixelU);
        printf("\nSaving I Color Buckets\n  ");
        pixelL.push_back(cb->min0);
        pixelU.push_back(cb->min0+CB0a-1);
        for (auto b : cb->bucket1) { save_bucket(b, rac, srcRanges, 1, pixelL, pixelU); pixelL[0] += CB0a; pixelU[0] += CB0a; }
        printf("\nSaving Q Color Buckets\n  ");
        pixelL[0] = cb->min0;
        pixelU[0] = cb->min0+CB0b-1;
        pixelL.push_back(cb->min1);
        pixelU.push_back(cb->min1+CB1-1);
        for (auto bv : cb->bucket2) {
                pixelL[1] = cb->min1;
                pixelU[1] = cb->min1+CB1-1;
                for (auto b : bv) {
                        save_bucket(b, rac, srcRanges, 2, pixelL, pixelU);
                        pixelL[1] += CB1; pixelU[1] += CB1;
                }
                pixelL[0] += CB0b; pixelU[0] += CB0b;
        }
        printf("\n");
        if (srcRanges->numPlanes() > 3) {
          printf("Saving Alpha Color Bucket: ");
          save_bucket(cb->bucket3, rac, srcRanges, 3, pixelL, pixelU);
          printf("\n");
        }
    }

    bool process(const ColorRanges *srcRanges, const Image &image) {
//            if (srcRanges->numPlanes() != 3) printf("OOPS, was expecting 3 planes!\n");
            std::vector<ColorVal> pixel;
            // fill buckets
            for (int r=0; r<image.rows(); r++) {
                for (int c=0; c<image.cols(); c++) {
                  pixel.clear();
                  for (int p=0; p<image.numPlanes(); p++) {
                    ColorVal v = image(p,r,c);
                    pixel.push_back(v);
                  }
                  cb->addColor(pixel);
                  if (totaldiscretecolors > 20000 && totalcontinuousbuckets > 2000) {
                        printf("Too many colors, not using color buckets.\n");
                        return false;
                  }
                }
            }

            cb->bucket0.simplify_lossless();
            cb->bucket3.simplify_lossless();


        // TODO: IMPROVE THESE HEURISTICS!
        // TAKE IMAGE SIZE INTO ACCOUNT!
        // CONSIDER RELATIVE AREA OF BUCKETS / BOUNDS!

            printf("Filled color buckets with %i discrete colors + %i continous buckets\n",totaldiscretecolors,totalcontinuousbuckets);

            if (totaldiscretecolors < 10000 && totalcontinuousbuckets < 500) return true;

            // simplify buckets
            for (auto& b : cb->bucket1) b.simplify(80);
            for (auto& bv : cb->bucket2) for (auto& b : bv) b.simplify(60);

            printf("Filled color buckets with %i discrete colors + %i continous buckets\n",totaldiscretecolors,totalcontinuousbuckets);
            if (totaldiscretecolors > 1000) {
                printf("Too many colors, simplifying...\n");
                for (auto& b : cb->bucket1) b.simplify(50);
                for (auto& bv : cb->bucket2) for (auto& b : bv) b.simplify(20);
                printf("Filled color buckets with %i discrete colors + %i continous buckets\n",totaldiscretecolors,totalcontinuousbuckets);
                if (totaldiscretecolors > 2000 || totalcontinuousbuckets > 2000) {
                  printf("Still too many colors, not using auto-indexing.\n");
                  return false;
                }
            }

            return true;
    }
};

#endif
