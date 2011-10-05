#include <string.h>
#include <stdint.h>

#include "chance.h"

void build_table(uint16_t *zero_state, uint16_t *one_state, size_t size, int factor, int max_p) {
  const int64_t one = 1LL << 32;
  int64_t p;
  int last_p8, p8, i;

  memset(zero_state,0,sizeof(uint16_t) * size);
  memset(one_state,0,sizeof(uint16_t) * size);

  last_p8 = 0;
  p = one / 2;
  for (i = 0; i < size / 2; i++) {
    p8 = (size * p + one / 2) >> 32; //FIXME try without the one
    if (p8 <= last_p8) p8 = last_p8 + 1;
    if (last_p8 && last_p8 < size && p8 <= max_p) one_state[last_p8] = p8;

    p += ((one - p) * factor + one / 2) >> 32;
    last_p8 = p8;
  }

  for (i = size - max_p; i <= max_p; i++) {
    if (one_state[i]) continue;

    p = (i * one + size / 2) / size;
    p += ((one - p) * factor + one / 2) >> 32;
    p8 = (size * p + one / 2) >> 32; //FIXME try without the one
    if (p8 <= i) p8 = i + 1;
    if (p8 > max_p) p8 = max_p;
    one_state[i] = p8;
  }

  for (i = 1; i < size; i++)
    zero_state[i] = size - one_state[size - i];

/*  for (int i=0; i<size; i++) {
    fprintf(stderr,"%i -> [%i,%i]\n",i,zero_state[i],one_state[i]);
  } */
}

SimpleBitChanceTable sbcTable(8);
