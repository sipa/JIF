#ifndef rac16_H_
#define rac16_H_ 1

#include <stdlib.h>
#include <stdint.h>

typedef uint32_t rac16_t;

typedef struct {
  rac16_t range;
  rac16_t low;
  int delayed_byte;
  int delayed_count;
  void *data;
  void (*writefn)(void *, int);
  int (*readfn)(void *);
} rac16_ctx_t;

void rac16_put_bit_frac(rac16_ctx_t *ctx, int num, int denom, int value);
void rac16_put_bit_b16(rac16_ctx_t *ctx, uint16_t b16, int value);

int rac16_get_bit_frac(rac16_ctx_t *ctx, int num, int denom);
int rac16_get_bit_b16(rac16_ctx_t *ctx, uint16_t b16);

void rac16_init_enc(rac16_ctx_t *ctx, void(*writefn)(void *, int), void *data);
void rac16_init_dec(rac16_ctx_t *ctx, int(*readfn)(void *), void *data);

void rac16_flush(rac16_ctx_t *ctx);

void rac16_build_table(uint16_t *zero_state, uint16_t *one_state, size_t size, int factor, int max_p);

#endif
