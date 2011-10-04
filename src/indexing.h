#ifndef _INDEXING_H_
#define _INDEXING_H_


#include "config.h"
#include "bitvector.h"
#include "image/pixel.h"

#define BUCKETS_Y 32
#define BUCKETS_I 7
#define BUCKETS_Q 7

#define MAX_PER_BUCKET 32

#define SIZE_Y (256/BUCKETS_Y)
#define SIZE_I (511/BUCKETS_I)
#define SIZE_Q (511/BUCKETS_Q)

typedef struct {
   uint16_t count;
   pixel_t color[MAX_PER_BUCKET];
} color_bucket;

typedef struct {
   color_bucket bucketYIQ[BUCKETS_Y][BUCKETS_I][BUCKETS_Q];
   int bucketYI[256][511];
   int bucketY[256];
   int bucketYIc[256][512];
   int bucketYc[257];
//   color_bucket bucket;
   int used;
   int max_diff[3];
} colors;

//int reduce_range(colors *a, int *min, int *max, int c, pixel_t *known, int qz, int guess);
//int reduce_range_Y(colors *a, int *y_min, int *y_max);
//int reduce_range_I(colors *a, int *i_min, int *i_max, int y);
//int reduce_range_Q(colors *a, int *q_min, int *q_max, int y, int i);
int same_color(pixel_t a, pixel_t b);
void round_color(const colors *a, pixel_t *pixel);

void add_color(colors *a, pixel_t *p);

void fill_bucketYI(colors *a, int glob_minq, int glob_maxq);


int bucket_exists(int yk, int ik, int qk);
void update_bounds(color_bucket *b, pixel_t *p);
void input_bounds(color_bucket *b, pixel_t *p, int n);
void init_bounds(color_bucket *b);
void init_bounds_input(color_bucket *b, int c);

//bitvector* compute_range(colors *a, int min, int max, int c, pixel_t *known, int qz, int guess, bitvector *v);
int compute_range(colors *a, const int min, const int max, const int c, const pixel_t *known, const int qz, const int guess, bitvector *bv);

void convert_indexed_clusters_to_full_buckets(colors *a);

int auto_indexing_heuristic(colors *a);


#endif
