#ifndef _CONTEXT_TYPES_H_
#define _CONTEXT_TYPES_H_


#include "config.h"
#include "symbol.h"
#include "color.h"
#include "indexing.h"
#include "image/pixel.h"
#include "image/image.h"


#define NB_PROPERTIES 11
#define NB_PROPERTIES_FFV1 10
#define NB_PROPERTIES_1 11
#define MAX_NB_PROPERTIES 11

#define NB_PROPERTIES_S 4


typedef struct {
   uint32_t count;
   int64_t sum[MAX_NB_PROPERTIES];
   symb_chs_t context;
   symb_chs_t splitcontext[MAX_NB_PROPERTIES*2];
   uint64_t splitbits[MAX_NB_PROPERTIES];
   uint64_t bits;
} context_leaf;

typedef struct {
   uint32_t count;
   int64_t sum[NB_PROPERTIES_S];
   chs_t context;
   chs_t splitcontext[NB_PROPERTIES_S*2];
   uint64_t splitbits[NB_PROPERTIES_S];
   uint64_t bits;
} context_leaf_bit;

typedef struct {
   int property; // 0 = leaf node
   union {
     context_leaf *context;
     struct {
       uint16_t branch;      // first child at this position, second child at branch+1
       uint16_t splitval;
       uint32_t count;
//       int64_t sum;
     } node; 
   } data;
} context_tree;

typedef struct {
   int property; // 0 = leaf node
   symb_chs_t *context;
   uint16_t branch;
   uint16_t splitval;
   uint32_t count;
   uint32_t nb;
   int64_t sum;
} context_tree_p2;


typedef struct {
   int property; // 0 = leaf node
   union {
     context_leaf_bit *context;
     struct {
       uint16_t branch;      // first child at this position, second child at branch+1
       uint16_t splitval;
//       uint32_t count;
//       int64_t sum;
     } node; 
   } data;
} context_tree_bit;

// full context table for range encoder
typedef struct {
// old static context structure:
//  chs_t splitContA[QZ_SIZE_S][QZ_DEPTH_S][3][4][QZ_BVAR_S]; // 0=split and continue      1=stop & interpolate
//  chs_t splitContL[QZ_SIZE_SL][QZ_DEPTH_SL][3][4][QZ_BVAR_S][QZ_LRDIFF_S];
//  symb_chs_t diff[QZ_SIZE][QZ_DEPTH][2][3][QZ_LRDIFF][QZ_GRAD][QZ_CROSS][QZ_VAR][QZ_BVAR][QZ_Y][QZ_I];

  context_tree_bit splitContA[3][MAX_CONTEXTS_S+2];
//  chs_t splitShort[32]; // 1=split along short axis  0=split along long axis and stop
//  chs_t splitAssym[32]; // 1=split assymetrically    0=split centrally
//  chs_t splitSided[32]; // 1=split near high coords  0=split near low coords
  chs_t interpolationMethod[QZ_SIZE_S][3][MAX_INTERPOLATIONS][MAX_INTERPOLATIONS]; 
  context_tree diff[3][MAX_CONTEXTS+2];
  context_tree_p2 diff_p2[3][MAX_CONTEXTS+2];
#if FULL_BORDERS
  context_tree_bit splitContL[3][MAX_CONTEXTS_S+2];
  context_tree diff1[3][MAX_CONTEXTS_1+2];
#endif  
  size_t treesize[4][3];                // [{diff/splitA/splitL/diff1}] [plane]
  symb_chs_t methods;
  symb_chs_t r_min[3];
  symb_chs_t r_max[3];

  symb_chs_t colors[3];
  symb_chs_t colors_Y;
  symb_chs_t colors_count0;
  symb_chs_t colors_count;

  symb_chs_t context_tree_property;
  symb_chs_t context_tree_count1;
  symb_chs_t context_tree_count2;

  uint64_t bitsSplit[4][2]; uint32_t symbSplit[4][2];
  uint64_t bitsInterpol[3]; uint32_t symbInterpol[3];
  uint64_t bitsPixeldata[3][2]; uint32_t symbPixeldata[3][2];
  uint64_t bitsContextTreeData[3];   
//  uint64_t bitsColorData;
  uint64_t bitsColorData[3];
  uint64_t bitsColorData_count;
  uint64_t bitsColorData_Y;
} ctx_t;

// data for encoder/decoder
typedef struct {
  image_t *source; // for encoder: source image
  image_t *rec; // image for reconstruction (during encoding and decoding)
//  uint64_t maxD[3]; // maximum color distances (encoder only)
  int maxD[3]; // maximum color distances (encoder only)
  int qz[3]; // quantization (1 for lossless)
  int sb[3]; // significant mantissa bits to be stored ({8,9,9} for lossless)
  double factor; // allowed factor (encoder only) area is OK if total dist < factor * maxD
  double factor_s; // small area factor (encoder only) per-pixel error is bounded by factor_s * maxD  (area of N pixels gets total dist < N*factor_s*maxD)
  symb_coder_t coder; // symbol encoder private data
  int range_min[3];  //
  int range_max[3];  // ranges for Y,I,Q
  int nb_pixels;

  colors color_data;

  ctx_t *ctx; // context table

  int traversal_method; // 0 = splitting,  1 = FFV1
  int phase; // 0 = all-in-one-go ; 1 = phase 1 ; 2 = phase 2

  int nb_interpolation_methods;
  int interpolation_methods[MAX_INTERPOLATIONS];

  // statistics
  int outputted_pixels[3]; // number of outputted pixels
  int splits[3][2]; // number of splits in each color plane
  int avg_Lsplit_depth[3];
  int avg_Asplit_depth[3];
  int pLeft,pRight,pTop,pBottom;
  int method[MAX_INTERPOLATIONS+1];
  uint64_t method_gain[MAX_INTERPOLATIONS+1];
  int minsplitA;
  int minsplitL;
} coder_t;




#endif
