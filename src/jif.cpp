
#include "maniac/rac.h"
#include "maniac/compound.h"
#include "maniac/util.h"

#include "image/image.h"

typedef std::vector<std::pair<int,int> > propRanges_t;
typedef std::vector<int> props_t;

void static initPropRanges(propRanges_t &propRanges, const Plane &plane)
{
    propRanges.clear();
    int min = plane.min - plane.max;
    int max = plane.max - plane.min;
    propRanges.push_back(std::make_pair(plane.min,plane.max));
    propRanges.push_back(std::make_pair(min,max));
    propRanges.push_back(std::make_pair(min,max));
    propRanges.push_back(std::make_pair(min,max));
    propRanges.push_back(std::make_pair(min,max));
    propRanges.push_back(std::make_pair(min,max));
/*
        properties[5] = quantize_log(UNSTRETCH(left->d[c])-UNSTRETCH(topleft->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[6] = quantize_log(UNSTRETCH(topleft->d[c])-UNSTRETCH(top->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[7] = quantize_log(UNSTRETCH(top->d[c])-UNSTRETCH(topright->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[8] = quantize_log(UNSTRETCH(toptop->d[c])-UNSTRETCH(top->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
        properties[9] = quantize_log(UNSTRETCH(leftleft->d[c])-UNSTRETCH(left->d[c]),(c==0 ? 255 : 510),SIGNED_LRDIFF,QZ_LRDIFF);
*/
}

void static calcProps(props_t &properties, const Plane &plane, int r, int c)
{
    properties.push_back(plane(r,c-1)-plane(r-1,c-1));  // left - topleft
    properties.push_back(plane(r-1,c-1)-plane(r-1,c));  // topleft - top
    properties.push_back(plane(r-1,c)-plane(r-1,c+1));  // top - topright
    properties.push_back(plane(r-2,c)-plane(r-1,c));    // toptop - top
    properties.push_back(plane(r,c-2)-plane(r,c-1));    // leftleft - left
}

template<typename I> void static swap(I& a, I& b)
{
    I c = a;
    a = b;
    b = c;
}

template<typename I> I static median3(I a, I b, I c)
{
    if (a<b) swap(a,b);
    if (b<c) swap(b,c);
    if (a<b) swap(a,b);
    return b;
}

ColorVal static predict(const Plane &plane, int r, int c)
{
    ColorVal left = plane(r-1,c);
    ColorVal top = plane(r,c-1);
    ColorVal gradient = left + top - plane(r-1,c-1);
    return median3(left,top,gradient);
}

typedef MultiscaleBitChance<SimpleBitChance> JifBitChance;

bool encode(const char* filename, const Image &image)
{
    FILE *f = fopen(filename,"w");
    RacOutput40 rac(f);

    SimpleSymbolCoder<JifBitChance, RacOutput40> metaCoder(rac, 24);
    int numPlanes = image.numPlanes();
    metaCoder.write_int(1, 16, numPlanes);
    for (int p = 0; p < numPlanes; p++) {
        const Plane& plane = image(p);
        metaCoder.write_int(1, 65536, plane.width);
        metaCoder.write_int(1, 65536, plane.height);
        metaCoder.write_int(-16777216, 16777215, plane.min);
        metaCoder.write_int(0, 16777216, plane.max - plane.min);
    }

    for (int p = 0; p < image.numPlanes(); p++) {
        const Plane &plane = image(p);
        printf("Plane #%i: %ix%i [%i..%i]\n", p, plane.width, plane.height, plane.min, plane.max);
        propRanges_t propRanges;
        initPropRanges(propRanges,plane);
        int nBits = ilog2((plane.max-plane.min)*2-1)+3;
        PropertySymbolCoder<JifBitChance, RacOutput40> coder(rac, propRanges, nBits);
        for (int r = 0; r < plane.height; r++) {
            for (int c = 0; c < plane.width; c++) {
                props_t properties;
                ColorVal guess = predict(plane,r,c);
                properties.push_back(guess);
                ColorVal curr = plane(r,c);
                calcProps(properties,plane,r,c);
                coder.write_int(properties, plane.min - guess, plane.max - guess, curr - guess);
//                fprintf(stderr, "%i(%i,%i)\n", curr - guess, plane.min - guess, plane.max - guess);
            }
        }
        printf("baldskjd\n");
    }

    rac.flush();
    fclose(f);
    return true;
}

bool decode(const char* filename, Image &image)
{
    image.reset();

    FILE *f = fopen(filename,"r");
    RacInput40 rac(f);

    SimpleSymbolCoder<JifBitChance, RacInput40> metaCoder(rac, 24);
    int numPlanes = metaCoder.read_int(1, 16);
    for (int p = 0; p < numPlanes; p++) {
        int width = metaCoder.read_int(1, 65536);
        int height = metaCoder.read_int(1, 65536);
        int min = metaCoder.read_int(-16777216, 16777215);
        int max = metaCoder.read_int(0, 16777216) + min;
        printf("Plane #%i: %ix%i [%i..%i]\n", p, width, height, min, max);
        image.add_plane(width, height, min, max);
    }

    for (int p = 0; p < image.numPlanes(); p++) {
        Plane &plane = image(p);
        propRanges_t propRanges;
        initPropRanges(propRanges,plane);
        int nBits = ilog2((plane.max-plane.min)*2-1)+3;
        PropertySymbolCoder<JifBitChance, RacInput40> coder(rac, propRanges, nBits);
        for (int r = 0; r < plane.height; r++) {
            for (int c = 0; c < plane.width; c++) {
                props_t properties;
                ColorVal guess = predict(plane,r,c);
                properties.push_back(guess);
                calcProps(properties,plane,r,c);
                ColorVal curr = coder.read_int(properties, plane.min - guess, plane.max - guess) + guess;
//                fprintf(stderr, "%i(%i,%i)\n", curr - guess, plane.min - guess, plane.max - guess);
                plane(r,c) = curr;
            }
        }
    }

    fclose(f);
    return true;
}


int main(int argc, char **argv)
{
    Image image;
    if (argc == 3) {
        image.load(argv[1]);
        encode(argv[2], image);
    } else if (argc == 4) {
        decode(argv[2], image);
        image.save(argv[3]);
    } else if (argc == 5) {
        image.load(argv[3]);
        image.save(argv[4]);
    } else {
        fprintf(stderr, "Usage: %s [-d] source dest\n", argv[0]);
        return 1;
    }
    return 0;
}
