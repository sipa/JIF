
#include "maniac/rac.h"
#include "maniac/compound.h"

#include "image/image.h"

typedef std::vector<std::pair<int,int> > propRanges_t;
typedef std::vector<int> props_t;
template<typename Rac> class coders_t : public std::vector<PropertySymbolCoder<SimpleBitChance, Rac> > {};

void static initPropRanges(propRanges_t &propRanges)
{
    propRanges.clear();
}

template<typename Rac> void static initCoders(const Image& image, Rac &rac, coders_t<Rac> &coders)
{
    for (int plane = 0; plane < image.numPlanes(); plane++)
    {
        propRanges_t propRanges;
        initPropRanges(propRanges);
        PropertySymbolCoder<SimpleBitChance, Rac> coder(rac,propRanges,8);
        coders.push_back(coder);
    }
}

bool encode(const char* filename, const Image &image)
{
    FILE *f = fopen(filename,"w");
    RacOutput40 rac(f);
    coders_t<RacOutput40> coders;
    initCoders(image, rac, coders);

    for (int p = 0; p < image.numPlanes(); p++) {
        const Plane &plane = image(p);
        ColorVal prev = 0;
        for (unsigned int r = 0; r < plane.height; r++) {
            for (unsigned int c = 0; c < plane.width; c++) {
                props_t properties;
                ColorVal curr = plane(r,c);
                coders[p].write_int(properties, -255, 255, default_range_test, curr - prev);
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
