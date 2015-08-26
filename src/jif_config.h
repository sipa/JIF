

// output the first K zoomlevels without building trees (too little data to learn much)
#define NB_NOLEARN_ZOOMS 12
// this is enough to get a reasonable thumbnail/icon before the tree gets built/transmitted

// more makes encoding more expensive, but results in better trees (smaller files)
#define TREE_LEARN_REPEATS 3


//#define PERMUTEPLANES 1

// if set, double the planes: most significant bits go first
//#define FIRSTQUANTIZE 1
#define QUANTIZATION 8
#define QUANTIZATIONSHIFT 3
#define QUANTIZATIONMASK 7
// drop the 3 least significant bits, e.g. rgb888 goes to rgb555 + rgb333, yiq899 goes to yiq566 + yiq333

#define CHECK_FOR_BROKENFILES 1

// broken:
//#define SMOOTHZOOM 1


#include "maniac/rac.h"
typedef RacInput40 RacIn;
typedef RacOutput40 RacOut;
//typedef RacInput24 RacIn;
//typedef RacOutput24 RacOut;
