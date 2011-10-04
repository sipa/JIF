#ifndef _SYMBOL_H_
#define _SYMBOL_H_ 1

#include <stdio.h>
#include "config.h"
#include "chance.h"
#include "bitvector.h"
//#include "indexing.h"

#if (RAC_BITS == 40)
#include "rac/rac40.h"
#define rac_ctx_t         rac40_ctx_t
#define rac_put_bit_b16   rac40_put_bit_b16
#define rac_get_bit_b16   rac40_get_bit_b16
#define rac_put_bit_frac  rac40_put_bit_frac
#define rac_get_bit_frac  rac40_get_bit_frac
#define rac_init_enc      rac40_init_enc
#define rac_init_dec      rac40_init_dec
#define rac_flush         rac40_flush
#endif

#if (RAC_BITS == 24)
#include "rac/rac24.h"
#define rac_ctx_t         rac24_ctx_t
#define rac_put_bit_b16   rac24_put_bit_b16
#define rac_get_bit_b16   rac24_get_bit_b16
#define rac_put_bit_frac  rac24_put_bit_frac
#define rac_get_bit_frac  rac24_get_bit_frac
#define rac_init_enc      rac24_init_enc
#define rac_init_dec      rac24_init_dec
#define rac_flush         rac24_flush
#endif

#if (RAC_BITS == 16)
#include "rac/rac16.h"
#define rac_ctx_t         rac16_ctx_t
#define rac_put_bit_b16   rac16_put_bit_b16
#define rac_get_bit_b16   rac16_get_bit_b16
#define rac_put_bit_frac  rac16_put_bit_frac
#define rac_get_bit_frac  rac16_get_bit_frac
#define rac_init_enc      rac16_init_enc
#define rac_init_dec      rac16_init_dec
#define rac_flush         rac16_flush
#endif

#if (ACCURATE_FRAC == 1)
// accurate symbol chance fractions use a separate counter for each possibility
typedef struct {
  int ch[1021];
}symb_chs_t;
#else
// normal symbol chance tables use exponent/mantissa/sign counters
typedef struct {
  chs_t chs[30];
} symb_chs_t;
#endif

typedef struct {
  chs_table_t table;
  rac_ctx_t rac;
} symb_coder_t;

void symb_chs_init(symb_chs_t *sc);
void symb_init_write(symb_coder_t *c, FILE* output, int cutoff);
void symb_init_read(symb_coder_t *c, FILE* input, int cutoff);
void symb_skip_bit(symb_coder_t *c, chs_t *chs, int val);
void symb_put_bit(symb_coder_t *c, chs_t *chs, int val, uint64_t *count);
void symb_put_bit_v(symb_coder_t *c, chs_t *chs, int val, uint64_t *count, int virtua);
int symb_get_bit_c(symb_coder_t *c, chs_t *chs, uint64_t *count);
int symb_get_bit(symb_coder_t *c, chs_t *chs);
void symb_put_simple_bit(symb_coder_t *c, int val, uint64_t *count);
int symb_get_simple_bit(symb_coder_t *c);

void symb_put_simple_int(symb_coder_t *c, int val, int min, int max, uint64_t *count);
int symb_get_simple_int(symb_coder_t *c, int min, int max);
void symb_put_simple_int_0(symb_coder_t *c, int val, int min, int max, uint64_t *count);
int symb_get_simple_int_0(symb_coder_t *c, int min, int max);

void symb_put_int(symb_coder_t *c, symb_chs_t *sc, int val, int min, int max, uint64_t *count);
int symb_get_int(symb_coder_t *c, symb_chs_t *sc, int min, int max);
void symb_put_big_int(symb_coder_t *c, symb_chs_t *sc1, symb_chs_t *sc2, int val, int min, int max, uint64_t *count);
int symb_get_big_int(symb_coder_t *c, symb_chs_t *sc1, symb_chs_t *sc2, int min, int max);
//void symb_put_int_limited_v_a(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, int virtua, bitvector *bv); //colors *a, int channel, pixel_t *pixel, int qz, int guess);
//void symb_put_int_limited_v_a(symb_coder_t *c, symb_chs_t *sc, int *val, const int min, const int max, uint64_t *count, const int nb_bits, const int virtua, bitvector *bv);
void symb_put_int_limited_a(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, bitvector *bv); //, colors *a, int channel, pixel_t *pixel, int qz, int guess);
void symb_put_int_limited_a1(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, bitvector *bv); //, colors *a, int channel, pixel_t *pixel, int qz, int guess);
void symb_put_int_limited_v(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, int virtua);
void symb_put_int_limited(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits);
int symb_get_int_limited_c_a(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits, uint64_t *count, bitvector *bv); //, colors *a, int channel, pixel_t *pixel, int qz, int guess);
int symb_get_int_limited_c(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits, uint64_t *count);
int symb_get_int_limited(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits);
void symb_flush(symb_coder_t *c);
void symb_show_stats(FILE *f, symb_coder_t *c);
void symb_cp(symb_chs_t *from,symb_chs_t *to);


#endif
