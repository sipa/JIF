#ifndef _BITVECTOR_H_
#define _BITVECTOR_H_ 1

#include <stdint.h>
#include <stdlib.h>
#include "assert.h"
#include "config.h"



typedef struct {
    int bits[1023];
//    int cumul[1024];
    int nb;
    int *pzero;
    int *pzeroc;
    int min;
    int max;
} bitvector;

#define EMPTY_VECTOR(v)  (( v==NULL || v->nb<0 ))





int get_first_in_range(const bitvector *v, int min, int max);

int get_bit(const bitvector *v, const int pos);

//void set_bit(bitvector *v, unsigned int pos);
//unsigned int count_range(bitvector *v, const int min, const int max, const int orig_min);
//unsigned int check_range(bitvector *v, const int min, const int max, const int orig_min);

//void disable_vector();
//void enable_vector();

static inline void init_vector(bitvector *v, int size) {
  assert(size<1023);
  v->nb = size-1;
  for (int i=0; i < size ; i++) {
      v->bits[i] = 0;
  }
}
static inline void init_vector_size(bitvector *v, int size) {
  assert(size<1023);
  v->nb = size-1;
// v->cumul[0]=0;
}


static inline unsigned int check_range(const bitvector *v, int min, int max) {
   if (EMPTY_VECTOR(v)) return (min>max? 0 : 1);
   if (max > v->max) max = v->max;
   if (min < v->min) min = v->min;
   if (max < min) return 0;
   if (v->pzeroc != NULL) return (v->pzeroc[max+1] - v->pzeroc[min]);
   for (int i=min ; i <= max ; i++) {
        if (v->pzero[i]) { return 1;}
   }

/*
   int begin=min-orig_min;
   int end=max-orig_min;
   if (begin<0) begin=0;
   if (end>v->nb) end=v->nb;
   if (end<begin) return 0;
//   return v->cumul[end+1] - v->cumul[begin];

   for (int i=begin ; i <= end ; i++) {
        if (v->bits[i]) { return 1;}
   }
   */
   return 0;
   
}

// no actual count needed, just 0, 1, or more (2)
static inline unsigned int count_range(const bitvector *v, int min, int max) {
   if (EMPTY_VECTOR(v)) return (min<max? 2 : (min==max? 1 : 0));
   if (max > v->max) max = v->max;
   if (min < v->min) min = v->min;
   if (max < min) return 0;
   if (v->pzeroc != NULL) return (v->pzeroc[max+1] - v->pzeroc[min]);
   unsigned int c = 0;
   for (int i=min ; i <= max ; i++) {
        if (v->pzero[i]) { 
         if (++c > 1) break;
        }
   }

/*
   int begin=min-orig_min;
   int end=max-orig_min;
   if (begin<0) begin=0;
   if (end>v->nb) end=v->nb;
   if (end<begin) return 0;
                // 0..end       0..begin-1
//   return v->cumul[end+1] - v->cumul[begin];

   unsigned int c = 0;
   for (unsigned int i = begin; i <= end ; i++) {
      if (v->bits[i]) {
         if (++c > 1) break;
      }
   }
   */
   return c;
   
}


static inline unsigned int check_position(const bitvector *v, const int p) {
//  return check_range(v,p,p,orig_min);

  if (EMPTY_VECTOR(v)) return 1;

  if (p > v->max) return 0;
  if (p < v->min) return 0;
  
  return (v->pzero[p] != 0);
  
/*  
  int pos=p-orig_min;
  if (pos < 0) return 0;
  if (pos > v->nb) return 0;
  return (v->bits[pos] != 0);
  */
  
}


static inline void set_bit(bitvector *v, unsigned int pos) {
  assert(pos <= v->nb);
  assert(pos >= 0);
  v->bits[pos] = 1;
//  assert(get_bit(v,pos));
}
static inline void change_bit(bitvector *v, unsigned int pos, int val) {
  assert(pos <= v->nb);
  assert(pos >= 0);
  v->bits[pos] = val;
//  v->cumul[pos+1] = v->cumul[pos]+val;
}

/*
static inline void compute_cumul(bitvector *v, const int min, const int max) {
  int count = 0;
  for (int i=min; i<=max; i++) {
        if (v->bits[i]) count++;
        v->cumul[i+1] = count;
  }
}
*/



#endif
