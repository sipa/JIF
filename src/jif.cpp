#include <string>

#include "maniac/rac.h"
#include "maniac/compound.h"
#include "maniac/util.h"

#include "image/image.h"
#include "image/color_range.h"
#include "transform/factory.h"

typedef std::vector<std::pair<int,int> > propRanges_t;
typedef std::vector<int> props_t;

void static initPropRanges(propRanges_t &propRanges, const ColorRanges &ranges, int p)
{
    propRanges.clear();
    int min = ranges.min(p);
    int max = ranges.max(p);
    int mind = min - max, maxd = max - min;

    propRanges.push_back(std::make_pair(min,max));

    for (int pp = 0; pp < p; pp++) {
        propRanges.push_back(std::make_pair(ranges.min(pp), ranges.max(pp)));
    }
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
    propRanges.push_back(std::make_pair(mind,maxd));
}


void static calcProps(props_t &properties, const Image &image, int p, int r, int c)
{
    for (int pp = 0; pp < p; pp++) {
        properties.push_back(image(pp,r,c));
    }
    properties.push_back(image(p,r,c-1)-image(p,r-1,c-1));  // left - topleft
    properties.push_back(image(p,r-1,c-1)-image(p,r-1,c));  // topleft - top
    properties.push_back(image(p,r-1,c)-image(p,r-1,c+1));  // top - topright
    properties.push_back(image(p,r-2,c)-image(p,r-1,c));    // toptop - top
    properties.push_back(image(p,r,c-2)-image(p,r,c-1));    // leftleft - left
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

typedef SimpleBitChance                         JifBitChancePass1;
typedef MultiscaleBitChance<3,SimpleBitChance>  JifBitChancePass2;
typedef MultiscaleBitChance<6,SimpleBitChance>  JifBitChanceTree;
typedef SimpleBitChance                         JifBitChanceMeta;

template<typename RAC> void static write_name(RAC& rac, std::string str)
{
    UniformSymbolCoder<RAC> coder(rac);
    coder.write_int(3, 8, str.size());
    for (unsigned int i=0; i<str.size(); i++) {
        char c = str[i];
        int n = ((c >= 'A' && c <= 'Z') ? c - 'A' :
                ((c >= 'a' && c <= 'z') ? c - 'a' :
                ((c >= '0' && c <= '9') ? c - '0' + 26 : 36)));
        coder.write_int(0, 36, n);
    }
}

template<typename RAC> std::string static read_name(RAC& rac)
{
    static char cs[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
    UniformSymbolCoder<RAC> coder(rac);
    int l = coder.read_int(3, 8);
    std::string str;
    for (int i=0; i<l; i++) {
        int n = coder.read_int(0, 36);
        str += cs[n];
    }
    return str;
}

ColorVal static predict(const Image &image, int p, int r, int c)
{
    ColorVal left = image(p,r-1,c);
    ColorVal top = image(p,r,c-1);
    ColorVal gradient = left + top - image(p,r-1,c-1);
    return median3(left,top,gradient);
}

template<typename Coder> void encode_ffv1_inner(std::vector<Coder*> &coders, const Image &image, const ColorRanges *ranges)
{
    for (int r = 0; r < image.rows(); r++) {
        for (int c = 0; c < image.cols(); c++) {
            for (int p = 0; p < image.numPlanes(); p++) {
                if (image.is_set(p,r,c)) {
                    props_t properties;
                    ColorVal guess = predict(image,p,r,c);
                    properties.push_back(guess);
                    ColorVal curr = image(p,r,c);
                    calcProps(properties,image,p,r,c);
                    coders[p]->write_int(properties, ranges->min(p,r,c) - guess, ranges->max(p,r,c) - guess, curr - guess);
                }
            }
        }
    }
}

template<typename Rac, typename Coder> void encode_ffv1_pass(Rac &rac, const Image &image, const ColorRanges *ranges, std::vector<Tree> &forest)
{
    std::vector<Coder*> coders;
    for (int p = 0; p < ranges->numPlanes(); p++) {
        propRanges_t propRanges;
        initPropRanges(propRanges, *ranges, p);
        int nBits = ilog2((ranges->max(p) - ranges->min(p))*2-1)+1;
        coders.push_back(new Coder(rac, propRanges, nBits, forest[p]));
    }

    encode_ffv1_inner(coders, image, ranges);

    for (int p = 0; p < image.numPlanes(); p++) {
        delete coders[p];
    }
}

template<typename BitChance, typename Rac> void encode_tree(Rac &rac, const ColorRanges *ranges, const std::vector<Tree> &forest)
{
    for (int p = 0; p < ranges->numPlanes(); p++) {
        propRanges_t propRanges;
        initPropRanges(propRanges, *ranges, p);
        MetaPropertySymbolCoder<BitChance, Rac> metacoder(rac, propRanges);
        metacoder.write_tree(forest[p]);
    }
}

bool encode(const char* filename, Image &image, std::vector<std::string> transDesc)
{
    FILE *f = fopen(filename,"w");
    RacOut rac(f);

    write_name(rac, "JIF1");

    SimpleSymbolCoder<JifBitChanceMeta, RacOut> metaCoder(rac, 24);
    int numPlanes = image.numPlanes();
    metaCoder.write_int(1, 16, numPlanes);
    metaCoder.write_int(1, 65536, image.cols());
    metaCoder.write_int(1, 65536, image.rows());
    for (int p = 0; p < numPlanes; p++) {
        metaCoder.write_int(1, 64, image.subSampleR(p));
        metaCoder.write_int(1, 64, image.subSampleC(p));
        metaCoder.write_int(-16777216, 16777215, image.min(p));
        metaCoder.write_int(0, 16777216, image.max(p) - image.min(p));
    }

    std::vector<const ColorRanges*> rangesList;
    std::vector<Transform*> transforms;
    rangesList.push_back(getRanges(image));
    for (unsigned int i=0; i<transDesc.size(); i++) {
        Transform *trans = create_transform(transDesc[i]);
        if (!trans->init(rangesList.back()) || !trans->process(rangesList.back(), image)) {
            fprintf(stderr, "Transform '%s' failed\n", transDesc[i].c_str());
        } else {
            printf("Doing transform '%s'\n", transDesc[i].c_str());
            rac.write(true);
            write_name(rac, transDesc[i]);
            trans->save(rangesList.back(), rac);
            rangesList.push_back(trans->meta(image, rangesList.back()));
            trans->data(image);
        }
    }
    rac.write(false);
    const ColorRanges* ranges = rangesList[1];

    // two passes
    std::vector<Tree> forest(ranges->numPlanes(), Tree());
    RacDummy dummy;
    fprintf(stdout,"Encoding data (pass 1)\n");
    encode_ffv1_pass<RacDummy, PropertySymbolCoder<JifBitChancePass1, RacDummy> >(dummy, image, ranges, forest);
    fprintf(stdout,"Encoding tree\n");
    encode_tree<JifBitChanceTree, RacOut>(rac, ranges, forest);
    fprintf(stdout,"Encoding data (pass 2)\n");
    encode_ffv1_pass<RacOut, FinalPropertySymbolCoder<JifBitChancePass2, RacOut> >(rac, image, ranges, forest);
    fprintf(stdout,"Encoding done\n");

    for (int i=transforms.size()-1; i>=0; i--) {
        delete transforms[i];
    }
    transforms.clear();
    for (unsigned int i=0; i<rangesList.size(); i++) {
        delete rangesList[i];
    }
    rangesList.clear();

    rac.flush();
    fclose(f);
    return true;
}

template<typename Coder> void decode_ffv1_inner(std::vector<Coder*> &coders, Image &image, const ColorRanges *ranges)
{
    for (int r = 0; r < image.rows(); r++) {
        for (int c = 0; c < image.cols(); c++) {
            for (int p = 0; p < image.numPlanes(); p++) {
                if (image.is_set(p,r,c)) {
                    props_t properties;
                    ColorVal guess = predict(image,p,r,c);
                    properties.push_back(guess);
                    calcProps(properties,image,p,r,c);
                    ColorVal curr = coders[p]->read_int(properties, ranges->min(p,r,c) - guess, ranges->max(p,r,c) - guess) + guess;
                    image(p,r,c) = curr;
                }
            }
        }
    }
}

template<typename BitChance, typename Rac> void decode_tree(Rac &rac, const ColorRanges *ranges, std::vector<Tree> &forest)
{
    for (int p = 0; p < ranges->numPlanes(); p++) {
        propRanges_t propRanges;
        initPropRanges(propRanges, *ranges, p);
        MetaPropertySymbolCoder<BitChance, Rac> metacoder(rac, propRanges);
        metacoder.read_tree(forest[p]);
    }
}

template<typename Rac, typename Coder> void decode_ffv1_pass(Rac &rac, Image &image, const ColorRanges *ranges, std::vector<Tree> &forest)
{
    std::vector<Coder*> coders;
    for (int p = 0; p < image.numPlanes(); p++) {
        propRanges_t propRanges;
        initPropRanges(propRanges, *ranges, p);
        int nBits = ilog2((ranges->max(p) - ranges->min(p))*2-1)+1;
        coders.push_back(new Coder(rac, propRanges, nBits, forest[p]));
    }

    decode_ffv1_inner(coders, image, ranges);

    for (int p = 0; p < image.numPlanes(); p++) {
        delete coders[p];
    }
}


bool decode(const char* filename, Image &image)
{
    image.reset();

    FILE *f = fopen(filename,"r");
    RacIn rac(f);

    std::string str = read_name(rac);
    if (str != "JIF1") {
        fprintf(stderr,"Unknown magic '%s'\n", str.c_str());
        return false;
    }

    SimpleSymbolCoder<JifBitChanceMeta, RacIn> metaCoder(rac, 24);
    int numPlanes = metaCoder.read_int(1, 16);
    int width = metaCoder.read_int(1, 65536);
    int height = metaCoder.read_int(1, 65536);
    image.init(width, height, 0, 0, 0);
    for (int p = 0; p < numPlanes; p++) {
        int subSampleR = metaCoder.read_int(1, 64);
        int subSampleC = metaCoder.read_int(1, 64);
        int min = metaCoder.read_int(-16777216, 16777215);
        int max = metaCoder.read_int(0, 16777216) + min;
        image.add_plane(min, max, subSampleR, subSampleC);
    }

    std::vector<const ColorRanges*> rangesList;
    std::vector<Transform*> transforms;
    rangesList.push_back(getRanges(image));
    while (rac.read()) {
        std::string desc = read_name(rac);
        Transform *trans = create_transform(desc);
        if (!trans) {
            fprintf(stderr,"Unknown transformation '%s'\n", desc.c_str());
            return false;
        }
        if (!trans->init(rangesList.back())) {
            fprintf(stderr,"Transformation '%s' failed\n", desc.c_str());
            return false;
        }
        printf("Doing transform '%s'\n", desc.c_str());
        trans->load(rangesList.back(), rac);
        rangesList.push_back(trans->meta(image, rangesList.back()));
        transforms.push_back(trans);
    }
    const ColorRanges* ranges = rangesList[1];

    std::vector<Tree> forest(ranges->numPlanes(), Tree());
    fprintf(stdout,"Decoding tree\n");
    decode_tree<JifBitChanceTree, RacIn>(rac, ranges, forest);
    fprintf(stdout,"Decoding data\n");
    decode_ffv1_pass<RacIn, FinalPropertySymbolCoder<JifBitChancePass2, RacIn> >(rac, image, ranges, forest);
    fprintf(stdout,"Decoding done\n");

    for (int i=transforms.size()-1; i>=0; i--) {
        transforms[i]->invData(image);
        delete transforms[i];
    }
    transforms.clear();
    for (unsigned int i=0; i<rangesList.size(); i++) {
        delete rangesList[i];
    }
    rangesList.clear();

    fclose(f);
    return true;
}


int main(int argc, char **argv)
{
    Image image;
    if (argc == 3) {
        image.load(argv[1]);
        std::vector<std::string> desc;
        desc.push_back("YIQ");
        desc.push_back("BND");
        encode(argv[2], image, desc);
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
