#ifndef _UTIL_H_
#define _UTIL_H_ 1

#include <stdint.h>

extern const uint8_t log2_tab[1024];
int ilog2(uint32_t l);
void indent(int n);

//#define STATS 1


#endif
