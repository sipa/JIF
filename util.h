#ifndef _UTIL_H_
#define _UTIL_H_ 1

#include <stdint.h>
#include "log4k.h"

extern const uint8_t log2_tab[1024];
int ilog2(uint32_t l);


int quantize_log(int val, int maxval, int sign, int levels);
int quantize_log_uint32(uint32_t val, int levels);
int qbvar(uint32_t var, int quant);

int static inline clamp(int val, int min, int max) {
  if (val<min) val=min;
  if (val>max) val=max;
  return val;
}

// division rounding down (instead of toward zero)
int static inline rdiv(int num, int den) {
  if (num>=0) {
    return num/den;
  } else {
    return (num-den+1)/den;
  }
}

#endif
