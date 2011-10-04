#ifndef _UTIL_H_
#define _UTIL_H_ 1

#include <stdint.h>
#include "log4k.h"
#include "assert.h"
#include <stdlib.h>


extern const uint8_t log2_tab[1024];
int ilog2(uint32_t l);



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

int static inline average(int64_t sum, uint32_t count) {
 assert(count>0);
// return (sum+count/2)/count;
 return sum/count;
}

int static inline quantize_log(int val, int maxval, int sign, int levels) {

/* faster version, not logarithmic

  assert(abs(val) <= maxval);
  if (sign == 0 && maxval<levels) return val;
//  if (sign == 1 && 2*maxval<levels) return val+maxval;
  
  if (levels==1 || val==0) return 0;
  int result;
  int rlevels=(sign ? (levels-1)/2 : levels-1);
  if (sign) {
    result = (abs(val)*(rlevels)/(maxval))*2 + (val<0);
  } else {
    result = abs(val)*(rlevels)/maxval+1;
  }
*/




  if (levels==1 || val==0) return 0;
  int rlevels=(sign ? (levels-1)/2 : levels-1);
  int pre=(log4k[abs(4096/val)]*rlevels)/(log4k[4096/maxval-1]);
  int result;  
  if (sign) {
    result = 1+pre*2+(val<0);  
  } else {
//    fprintf(stderr,"val:%i, maxval:%i, levels:%i, pre:%i\n",val,maxval,levels,pre);
    result = 1+pre;  
  }

  
  
  if (result >= levels) result=levels-1;
  assert(result < levels);
  assert(result >= 0);
//  if (result < 0 || result >= levels) fprintf(stderr,"Warning: quantized %i to %i which is not in target range 0..%i\n",val,result,levels);
  return result;
}


int median7(int a, int b, int c, int d, int e, int f, int g);
int median5(int a, int b, int c, int d, int e);
int median_3(int a, int b, int c);

#endif
