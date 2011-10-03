#ifndef RAC24_H_
#define RAC24_H_ 1

#include <stdlib.h>
#include <stdint.h>

typedef uint32_t rac24_t;

typedef struct {
  rac24_t range;
  rac24_t low;
  int delayed_byte;
  int delayed_count;
  void *data;
  void (*writefn)(void *, int);
  int (*readfn)(void *);
} rac24_ctx_t;

void rac24_put_bit_frac(rac24_ctx_t *ctx, int num, int denom, int value);
void rac24_put_bit_b16(rac24_ctx_t *ctx, uint16_t b16, int value);

int rac24_get_bit_frac(rac24_ctx_t *ctx, int num, int denom);
int rac24_get_bit_b16(rac24_ctx_t *ctx, uint16_t b16);

void rac24_init_enc(rac24_ctx_t *ctx, void(*writefn)(void *, int), void *data);
void rac24_init_dec(rac24_ctx_t *ctx, int(*readfn)(void *), void *data);

void rac24_flush(rac24_ctx_t *ctx);

void rac24_build_table(uint16_t *zero_state, uint16_t *one_state, size_t size, int factor, int max_p);

#endif
