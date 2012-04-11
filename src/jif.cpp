
#include "maniac/rac.h"
#include "maniac/compound.h"
#include "maniac/util.h"

#include "image/image.h"

typedef std::vector<std::pair<int,int> > propRanges_t;
typedef std::vector<int> props_t;

void static initPropRanges(propRanges_t &propRanges)
{
    propRanges.clear();
}

bool encode(const char* filename, const Image &image)
{
    FILE *f = fopen(filename,"w");
    RacOutput40 rac(f);

    SimpleSymbolCoder<SimpleBitChance, RacOutput40> metaCoder(rac, 24);
    int numPlanes = image.numPlanes();
    writer(metaCoder, 1, 16, default_range_test, numPlanes);
    for (int p = 0; p < image.numPlanes(); p++) {
        const Plane& plane = image(p);
        writer(metaCoder, 1, 65536, default_range_test, plane.width);
        writer(metaCoder, 1, 65536, default_range_test, plane.height);
        writer(metaCoder, -16777216, 16777215, default_range_test, plane.min);
        writer(metaCoder, 0, 16777216, default_range_test, plane.max - plane.min);
    }

    for (int p = 0; p < image.numPlanes(); p++) {
        const Plane &plane = image(p);
        propRanges_t propRanges;
        initPropRanges(propRanges);
        int nBits = ilog2((plane.max-plane.min)*2-1)+3;
        PropertySymbolCoder<SimpleBitChance, RacOutput40> coder(rac, propRanges, nBits);
        ColorVal prev = 0;
        for (unsigned int r = 0; r < plane.height; r++) {
            for (unsigned int c = 0; c < plane.width; c++) {
                props_t properties;
                ColorVal guess = prev;
                ColorVal curr = plane(r,c);
                coder.write_int(properties, plane.min - guess, plane.max - guess, default_range_test, curr - guess);
                prev = curr;
            }
        }
    }

    rac.flush();
    fclose(f);
    return true;
}

bool decode(const char* filename, Image &image)
{
    return true;
}


int main(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Usage: %s image.ppm image.jif\n", argv[0]);
        return 1;
    }
    Image image;
    image.load(argv[1]);
    encode(argv[2], image);
    return 0;
}
