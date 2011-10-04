#ifndef _CONTEXTS_H_
#define _CONTEXTS_H_

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <stdint.h>
//#include <math.h>

#include "bitvector.h"
//#include "config.h"
//#include "util.h"
//#include "symbol.h"
//#include "image/pixel.h"
//#include "image/image.h"
//#include "color.h"
//#include "interpol.h"
//#include "crc32k.h"
//#include "log4k.h"
//#include "indexing.h"
#include "context_types.h"

//#define DEBUGMODE




//#include "init_chances.h"

void init_context_leaf_counts(context_tree *first);
void init_context_leaf_splitcontexts(context_tree *first);
void init_context_leaf(context_tree *first);
void init_context_leaf_counts_bit(context_tree_bit *first);
void init_context_leaf_splitcontexts_bit(context_tree_bit *first);
void init_context_leaf_bit(context_tree_bit *first);
void ctx_init(ctx_t *ctx, int cutoff, int phase);

void output_context_tree(context_tree *tree, coder_t *coder, int treenumber, int c);
void output_context_tree_p2(context_tree_p2 *tree, coder_t *coder, int treenumber, int c);

void input_context_tree_p2(context_tree_p2 *tree, coder_t *coder, int treenumber, int c);
context_leaf* find_context(context_tree *tree,int *properties,coder_t *coder,int treenumber, int c, int *vcontext);
context_leaf_bit* find_context_bit(context_tree_bit *tree,int properties[NB_PROPERTIES_S],coder_t *coder,int treenumber, int c, int *vcontext);


void output_interpolation_methods(coder_t *encode);
void input_interpolation_methods(coder_t *decode);
void output_interpolation_choice(coder_t *encode, int i, int c, int method, int prev);
int input_interpolation_choice(coder_t *decode, int i, int c, int prev);
void update_avg_split_depths(coder_t *encode, int c, int depth, int line);

void output_mask(coder_t *encode,int oldmask,int newmask,int qsize,int depth,uint32_t *bvar, uint32_t *bdiffE, int line, const pixel_t *left, const pixel_t *right);
int input_mask(coder_t *decode,int mask,int qsize,int depth,uint32_t *bvar, uint32_t *bdiffE, int line, const pixel_t *left, const pixel_t *right);
void output_min_max(coder_t *encode,int min,int max,int oldmin,int oldmax, int c);
void input_min_max(coder_t *decode,int *min,int *max,int oldmin,int oldmax, int c);
void output_pixel(coder_t *encode, int x, int y, const pixel_t *guess, uint32_t size, int mask, 
                                const pixel_t *topleft, const pixel_t *topright, const pixel_t *left,const pixel_t *right, uint32_t bvar[3], uint32_t bdiffE[3], int depth, int onepixel, int splitdir);
void input_pixel(coder_t *decode, int x, int y, const pixel_t *guess, uint32_t size, int mask,
                                const pixel_t *topleft, const pixel_t *topright, const pixel_t *left,const pixel_t *right, uint32_t bvar[3], uint32_t bdiffE[3], int depth, int onepixel, int splitdir);

#if FFV1
void output_pixel_ffv1(coder_t *encode, int x, int y, const pixel_t *guess, int mask,
                                const pixel_t *topleft, const pixel_t *left, const pixel_t *top,const pixel_t *leftleft,const pixel_t *toptop,const pixel_t *topright);
void input_pixel_ffv1(coder_t *decode, int x, int y, const pixel_t *guess, int mask,
                                const pixel_t *topleft, const pixel_t *left, const pixel_t *top,const pixel_t *leftleft,const pixel_t *toptop,const pixel_t *topright);
#endif

void output_colors(coder_t *coder);
void input_colors(coder_t *coder);


#endif
