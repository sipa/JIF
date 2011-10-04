#ifndef _CHANCE_H_
#define _CHANCE_H_ 1

#include <stdint.h>
#include <stdio.h>

#include "config.h"
#include "log4k.h"

typedef struct {
#if (SYMBOL_STATS == 1)
   uint32_t count;
   uint32_t total;
#endif
#if (MULTISCALE_CHANCES == 1)
   uint32_t qual[MULTISCALE_LEVELS]; // estimate of number of bits produced, if corresponding average was used
   uint16_t chance[MULTISCALE_LEVELS]; // chance measured by each average
   uint16_t best;
#else
   uint16_t chance; // use 16-bit integers to represent chances (only 12 bits used)
#endif
} chs_t;

typedef struct {
  int cutoff;
#if (MULTISCALE_CHANCES == 1)
//  uint32_t count[MULTISCALE_LEVELS]; // how often each scale was used
  uint16_t next[MULTISCALE_LEVELS][2][4096]; // for each scale, a table to compute the next chance for 0 and 1
#else
  uint16_t next[2][4096]; // just two tables
#endif
} chs_table_t;

void static inline chs_init(chs_t *chs, int chance) {
#if (MULTISCALE_CHANCES == 1)
  for (int i=0; i<MULTISCALE_LEVELS; i++) chs->chance[i]=chance;
  for (int i=0; i<MULTISCALE_LEVELS; i++) chs->qual[i]=0;
  chs->best=0;
#else
  chs->chance=chance;
#if (SYMBOL_STATS == 1)
  chs->count=0;
  chs->total=0;
#endif
#endif
}

void static inline chs_cp(chs_t *from, chs_t *to) {
#if (MULTISCALE_CHANCES == 1)
  for (int i=0; i<MULTISCALE_LEVELS; i++) to->chance[i]=from->chance[i];
  for (int i=0; i<MULTISCALE_LEVELS; i++) to->qual[i]=from->qual[i];
  to->best = from->best;
#else
  to->chance=from->chance;
#if (SYMBOL_STATS == 1)
  to->count=from->count;
  to->total=from->total;
#endif
#endif
}

int static inline chs_get(chs_t *chs, chs_table_t *tbl) {
#if (MULTISCALE_CHANCES == 1)
/*  int best=0;
  for (int i=1; i<MULTISCALE_LEVELS; i++) { // try all levels
    if (chs->qual[i]<chs->qual[best]) best=i; // find the one with the best bits estimate
  }*/
  int best=chs->best;
//  assert(best>=0 && best <MULTISCALE_LEVELS);
//  tbl->count[best]++; // update statistics in table
  return chs->chance[best]; // return chance for best one
#else
#if (SYMBOL_STATS == 1)
#if (USE_AVERAGE_CHANCES == 1)
  if (chs->total>240) {
    uint64_t cnt=chs->count;
    uint64_t tot=chs->total;
    int ret= (((cnt<<13)+tot)/(tot*2));
    if (ret>4096-tbl->cutoff) ret=4096-tbl->cutoff;
    if (ret<tbl->cutoff) ret=tbl->cutoff;
    return ret;
  }
#endif
#endif
  return chs->chance;
#endif
}

void static inline chs_put_scales(chs_t *chs, chs_table_t *tbl, int val, int levels) {
#if (MULTISCALE_CHANCES == 1)
//  chs->best=0;
  for (unsigned int i=0; i<levels; i++) { // for each scale
    uint64_t sbits=log4k[val ? chs->chance[i] : 4096-chs->chance[i]]; // number of bits if this scale was used
    uint64_t oqual=chs->qual[i]; // previous estimate of bits used for this scale
//    chs->qual[i]=(oqual*4095+sbits*65537+2048)/4096; // update bits estimate (([0-2**32-1]*4095+[0..2**16-1]*65537+2048)/4096 = [0..2**32-1])
    chs->qual[i]=(oqual*4095+sbits*65537+2048)>>12; // update bits estimate (([0-2**32-1]*4095+[0..2**16-1]*65537+2048)/4096 = [0..2**32-1])
    if (chs->qual[i] < chs->qual[chs->best]) chs->best=i;
    chs->chance[i] = tbl->next[i][val][chs->chance[i]]; // update chance
  }
#else
  chs->chance = tbl->next[val][chs->chance];
#endif
#if (SYMBOL_STATS == 1)
  chs->count+=val;
  chs->total++;
#endif
}
void static inline chs_put(chs_t *chs, chs_table_t *tbl, int val) {
  chs_put_scales(chs, tbl, val, MULTISCALE_LEVELS);
}

void chs_table_init(chs_table_t *tbl, int offset);

void chs_show_stats(FILE *f, chs_table_t *tbl);

#endif
