#ifndef CRC32K_H_
#define CRC32K_H_ 1

#include <stdint.h>

extern uint_fast32_t const crc32k_tab[256];

#define crc32k_transform(var,byte) do { (var) = ((var) << 8) ^ crc32k_tab[(((var) >> 24) ^ (byte)) & 0xFF]; } while(0);

#endif
