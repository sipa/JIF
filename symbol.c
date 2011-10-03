#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "config.h"
#include "util.h"
#include "symbol.h"
#include "log4k.h"
#include "chance.h"

// output a bit to the range encoder, updating bit chance tables
void inline symb_put_bit(symb_coder_t *c, chs_t *chs, int val, uint64_t *count) {
  int chance=chs_get(chs,&c->table);
  rac_put_bit_b16(&c->rac,chance << 4,val);
  if (count) (*count) += log4k[val ? chance : 4096-chance];
  chs_put(chs,&c->table,val);
//  fprintf(stderr,"symb_put_bit: %i with chance %i\n",val,chance);
}

// input a bit from the range decoder, updating bit chance tables
int inline symb_get_bit(symb_coder_t *c, chs_t *chs) {
  int chance=chs_get(chs,&c->table);
  int val = rac_get_bit_b16(&c->rac,chance << 4);
  chs_put(chs,&c->table,val);
//  fprintf(stderr,"symb_get_bit: %i with chance %i\n",val,chance);
  return val;
}

void inline symb_put_simple_bit(symb_coder_t *c, int val, uint64_t *count) {
  rac_put_bit_b16(&c->rac,0x8000,val);
  if (count) (*count)++;
//  fprintf(stderr,"symb_put_bit: %i with chance simple\n",val);
}

int inline symb_get_simple_bit(symb_coder_t *c) {
  int val=rac_get_bit_b16(&c->rac,0x8000);
//  fprintf(stderr,"symb_get_bit: %i with chance simple\n",val);
  return val;
}

// do not output/input anything, but update chance tables anyway
void inline symb_skip_bit(symb_coder_t *c, chs_t *chs, int val) {
  chs_put(chs,&c->table,val);
}

#if (ACCURATE_FRAC == 1)
// initializer of symbol chance tables when using accurate symbol chance fractions
void symb_chs_init(symb_chs_t *sc) {
  for (int i=0; i<1021; i++) {
    int diff=abs(i-510);
    sc->ch[i]=10 + (1 << (8-log2_tab[diff]));
  }
}
#else
// initializer of symbol chance tables when using exponent/mantissa/sign representation
void symb_chs_init(symb_chs_t *sc) {
  for (int i = 0; i < 30; i++) {
    chs_init(&sc->chs[i],0x800);
  }
  
  // is zero?
  chs_init(&sc->chs[0],1500);

  // exponent
  chs_init(&sc->chs[1],3200);
  chs_init(&sc->chs[2],2800);
  chs_init(&sc->chs[3],2600);
  chs_init(&sc->chs[4],2400);
  chs_init(&sc->chs[5],2000);
  chs_init(&sc->chs[6],1500);
  chs_init(&sc->chs[7],800);
  chs_init(&sc->chs[8],300);
  
  // mantissa
  chs_init(&sc->chs[21],1800);
  chs_init(&sc->chs[22],1800);
  chs_init(&sc->chs[23],1800);
  chs_init(&sc->chs[24],1700);
  chs_init(&sc->chs[25],1600);
  chs_init(&sc->chs[26],1200);
  chs_init(&sc->chs[27],1000);
  chs_init(&sc->chs[28],800);
 
}
#endif

static int nbytes = 0;

int static cb_read(void* input) {
  int r = fgetc((FILE*) input);
  nbytes++;
  if (r < 0) return 0;
  return r;
}

void static cb_write(void* output, int data) {
  nbytes++;
  fputc(data,(FILE*) output);
}

void symb_init_read(symb_coder_t *c, FILE* input, int cutoff) {
  chs_table_init(&c->table,cutoff);
  rac_init_dec(&c->rac,cb_read,input);
  nbytes=0;
}

void symb_init_write(symb_coder_t *c, FILE* output, int cutoff) {
  chs_table_init(&c->table,cutoff);
  rac_init_enc(&c->rac,cb_write,output);
  nbytes=0;
}

// statistics
static uint64_t pos = 0, neg = 0, neut = 0;
static int64_t nb_diffs = 0, tot_diffs = 0, tot_sq_diffs = 0;

#if (ACCURATE_FRAC == 1)
// write out a single bit to the range encoder when using accurate fractions
void static put_sbit(symb_coder_t *c, int *chs, int num, int val, int sum, uint64_t *count) {
  if (num<=1) return;
  int med=num/2;
  int hsum=0;
  for (int i=med; i<num; i++) {
    hsum += chs[i];
  }
  //  fprintf(stderr,"put_sbit: %i/%i: %i (val=%i med=%i)\n",hsum,sum,val>=med,val,med);
  if (val>=med) {
    rac_put_bit_frac(&c->rac,hsum,sum,1);
    if (count) (*count) += log((double)hsum/(double)sum)*(-1.442695040888963407359924681);
    put_sbit(c,chs+med,num-med,val-med,hsum,count);
  } else {
    rac_put_bit_frac(&c->rac,hsum,sum,0);
    if (count) (*count) += log((double)lsum/(double)sum)*(-1.442695040888963407359924681);
    int lsum=sum-hsum;
    put_sbit(c,chs,med,val,lsum,count);
  }
}

// read in a single bit from the range encodering when using accurate fractions
int static get_sbit(symb_coder_t* c, int *chs, int num, int sum, int add) {
  if (num<=1) return add;
  int med=num/2;
  int hsum=0;
  for (int i=med; i<num; i++) {
    hsum += chs[i];
  }
  int bit=rac_get_bit_frac(&c->rac,hsum,sum);
  //  fprintf(stderr,"get_sbit: %i/%i: %i (med=%i)\n",hsum,sum,bit,med);
  if (bit) {
    return get_sbit(c,chs+med,num-med,hsum,add+med);
  } else {
    int lsum=sum-hsum;
    return get_sbit(c,chs,med,lsum,add);
  }
}

// update accurate fraction chance tables
void static update_chs(symb_chs_t *sc, int val) {
  for (int i=0; i<1021; i++) {
    sc->ch[i] = (sc->ch[i]*19+10)/20;
  }
  sc->ch[val] += 1000000;
}

// write out a whole symbol to the range encoder using accurate fractions
void symb_put_int(symb_coder_t *c, symb_chs_t *sc, int val, int min, int max) {
  int sum=0;
  assert(max>=min);
  assert(val>=min);
  assert(val<=max);
  for (int i=min; i<=max; i++) {
    sum += sc->ch[i+510];
  }
  put_sbit(c,sc->ch+min+510,max-min+1,val-min,sum);
  update_chs(sc,val+510);
}

// read in a whole symbol from the range encoder using accurate fractions
int symb_get_int(symb_coder_t *c, symb_chs_t *sc, int min, int max) {
  int sum=0;
  assert(max>=min);
  for (int i=min; i<=max; i++) {
    sum += sc->ch[i+510];
  }
  int val=get_sbit(c,sc->ch+min+510,max-min+1,sum,min);
  assert(val>=min);
  assert(val<=max);
  update_chs(sc,val+510);
  return val;
}

#else

/* write out a whole symbol to the range encoder using exponent/mantissa/sign representation
 - a zero bit is written if zero is within [min..max]
 - exponent bits are written in unary notation as long as no overflow is implied by them
 - mantissa bits are written only when necessary
 - a sign bit is written if both val and -val are within [min..max]
 */


inline void symb_put_int_limited(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits) {
//  fprintf(stderr,"jif_put_symb: %4i in [%4i..%4i]\n",val,min,max);
  // initialization
  assert(min<=max);
  assert(*val>=min);
  assert(*val<=max);
  if (min == max) return;

  int amax = abs(max) > abs(min) ? abs(max) : abs(min);
  int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));

  int mymax = 0;
  if (min*max < 0) {
        mymax = (*val>0 ? abs(max) : abs(min));
  } else {
        mymax = amax;
  }

  
  int bmax = log2_tab[amax];
  int bmin = log2_tab[amin];

  // statistics
  if (*val == 0) {
    neut++;
  } else if (*val > 0) {
    pos++;
  } else {
    neg++;
  }
  nb_diffs++;
  tot_diffs += *val;
  tot_sq_diffs += *val * *val;

  if (*val) { // nonzero
    const int a = abs(*val);
    const int e = log2_tab[a];

    // zeroness (chance state 0)
    if (amin == 0) { // zero is possible
      symb_put_bit(c,&sc->chs[0],0,count);
    }

    // unary encoding of exponent (chance states 1..9)
    assert(e<=9);
    int i = bmin;
    while (i < e) {
      symb_put_bit(c,&sc->chs[1 + i],1,count);
      i++;
    }
    if (e < bmax) symb_put_bit(c,&sc->chs[1 + i],0,count);

    // sign (chance states 11..20)
    if (min*max<0) {
      symb_put_bit(c,&sc->chs[11 + e],*val < 0,count);
    }


    // mantissa (chance states 21..29)
    int run = (1 << e);
    int left = (1 << e) - 1;
    for (int i = e-1; i>=0; i--) {
      left ^= (1 << i);
      int bit = 1;
      if (run + (1 << i) > mymax) { // 1-bit would cause overflow
        bit = 0;
        symb_skip_bit(c,&sc->chs[21 + i],0);
      } else if (run + left < amin) { // 0-bit would cause underflow
        symb_skip_bit(c,&sc->chs[21 + i],1);
      } else { // both 0 and 1 are possible
        if (i >= e-nb_bits) {
          bit = (a >> i) & 1;
          symb_put_bit(c,&sc->chs[21 + i],bit,count);
        } else if (i == e-nb_bits-1) {
          bit = 1;
        } else {
          bit = 0;
        }
      }
      run |= (bit << i);
    }
    *val = (*val<0?-run:run);


  } else { // zero
    symb_put_bit(c,&sc->chs[0],1,count);
  }
  assert(*val>=min);
  assert(*val<=max);
}
void symb_put_int(symb_coder_t *c, symb_chs_t *sc, int val, int min, int max, uint64_t *count) {
   symb_put_int_limited(c,sc, &val, min, max, count, 64);
}

// read in a whole symbol from the range encoder using exponent/mantissa/sign representation

inline int symb_get_int_limited(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits) {
  // initialization
  assert(min<=max);
  if (min == max) return min;

  int amax = abs(max) > abs(min) ? abs(max) : abs(min);
  int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
  int bmax = log2_tab[amax];
  int bmin = log2_tab[amin];

  int mymax = 0;

  if (amin == 0 && symb_get_bit(c,&sc->chs[0])) { // zero
//    fprintf(stderr,"jif_get_symb:    0 in [%4i..%4i]\n",min,max);
    return 0;
  } else { // nonzero
    // unary encoding of exponent (chance states 1..10)
    int e = bmin;
    while (e < bmax && symb_get_bit(c,&sc->chs[1 + e]))
      e++;
    assert(e<=9);

    // sign (chance states 11..20)
    int sign = 0;
    if (min*max < 0) {
      if (symb_get_bit(c,&sc->chs[11 + e])) {
        sign = -1;
      } else {
        sign = 1;
      }
    } else {
      sign = 1;
      if (max<=0) sign = -1;
    }


    if (min*max < 0) {
          mymax = (sign>0 ? abs(max) : abs(min));
    } else {
          mymax = amax;
    }

    // mantissa (chance states 21..29)
    int run = (1 << e);
    int left = (1 << e) - 1;
    for (int i = e-1; i>=0; i--) {
      left ^= (1 << i);
      int bit = 1;
      if (run + (1 << i) > mymax) { // 1-bit would cause overflow
        bit = 0;
        symb_skip_bit(c,&sc->chs[21 + i],0);
      } else if (run + left < amin) { // 0-bit would cause underflow
        symb_skip_bit(c,&sc->chs[21 + i],1);
      } else { // both 0 and 1 are possible
        if (i >= e-nb_bits) {
          bit = symb_get_bit(c,&sc->chs[21 + i]);
        } else if (i == e-nb_bits-1) {
          bit = 1;
        } else {
          bit = 0;
        }
      }
      run |= (bit << i);
    }
    int ret = run*sign;

    // output
    if (ret<min || ret>max) fprintf(stderr,"jif_get_symb: %4i in [%4i..%4i]\n",ret,min,max);
    assert(ret>=min);
    assert(ret<=max);
    return ret;
  }
}
int symb_get_int(symb_coder_t *c, symb_chs_t *sc, int min, int max) {
 return symb_get_int_limited(c, sc, min, max, 64);
}

#endif

void symb_flush(symb_coder_t *c) {
  rac_flush(&c->rac);
}

void symb_show_stats(FILE *f, symb_coder_t *c) {
/*  fprintf(f,"positive diffs: %llu\n",(unsigned long long) pos);
  fprintf(f,"negative diffs: %llu\n",(unsigned long long) neg);
  fprintf(f,"neutral  diffs: %llu\n",(unsigned long long) neut);
  double avg = 1.0 * tot_diffs / nb_diffs;
  double dev = sqrt(1.0 * tot_sq_diffs / nb_diffs - avg * avg);
  fprintf(f,"diff total: %lli / nb_diffs: %llu  = avg diff: %g (+-%g)\n",(long long) tot_diffs,
      (unsigned long long) nb_diffs,avg,dev); */
  fprintf(f,"bytes read/written: %i\n",nbytes);
  chs_show_stats(f,&c->table);
}
