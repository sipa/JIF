#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "config.h"
#include "util.h"
#include "symbol.h"
#include "log4k.h"
#include "chance.h"
//#include "indexing.h"
#include "bitvector.h"

// output a bit to the range encoder, updating bit chance tables
void inline symb_put_bit(symb_coder_t *c, chs_t *chs, int val, uint64_t *count) {
  int chance=chs_get(chs,&c->table);
  rac_put_bit_b16(&c->rac,chance << 4,val);
  if (count) (*count) += log4k[val ? chance : 4096-chance];
  chs_put(chs,&c->table,val);
//  fprintf(stderr,"symb_put_bit: %i with chance %i\n",val,chance);
}

// output a bit to the range encoder, updating bit chance tables
void inline symb_put_bit_v(symb_coder_t *c, chs_t *chs, int val, uint64_t *count, int virtua) {
  int chance=chs_get(chs,&c->table);
  if (!virtua) rac_put_bit_b16(&c->rac,chance << 4,val);
  if (count) (*count) += log4k[val ? chance : 4096-chance];
  chs_put(chs,&c->table,val);
//  fprintf(stderr,"symb_put_bit: %i with chance %i\n",val,chance);
}

// input a bit from the range decoder, updating bit chance tables
int inline symb_get_bit_c(symb_coder_t *c, chs_t *chs, uint64_t *count) {
  int chance=chs_get(chs,&c->table);
  int val = rac_get_bit_b16(&c->rac,chance << 4);
  if (count) (*count) += log4k[val ? chance : 4096-chance];
  chs_put(chs,&c->table,val);
//  fprintf(stderr,"symb_get_bit: %i with chance %i\n",val,chance);
  return val;
}
int inline symb_get_bit(symb_coder_t *c, chs_t *chs) {
 return symb_get_bit_c(c, chs,NULL);
}
void inline symb_put_simple_bit(symb_coder_t *c, int val, uint64_t *count) {
  rac_put_bit_b16(&c->rac,0x8000,val);
  if (count) (*count) += log4k[2048];
//  if (count) (*count)++;
//  fprintf(stderr,"symb_put_bit: %i with chance simple\n",val);
}

int inline symb_get_simple_bit(symb_coder_t *c) {
  int val=rac_get_bit_b16(&c->rac,0x8000);
//  fprintf(stderr,"symb_get_bit: %i with chance simple\n",val);
  return val;
}


void inline symb_put_simple_int(symb_coder_t *c, int val, int min, int max, uint64_t *count) {
  assert(min <= max);
  assert(val >= min);
  assert(val <= max);
  int x = val-min;
  int bits = ilog2(max-min)+1;
  for (int i=0; i<bits; i++) {
        symb_put_simple_bit(c, (x >> i) & 1, count);
  }
}
int inline symb_get_simple_int(symb_coder_t *c, int min, int max) {
  assert(min <= max);
  int x = 0;
  int bits = ilog2(max-min)+1;
  for (int i=0; i<bits; i++) {
        x |= (symb_get_simple_bit(c) << i);
  }
  return x+min;
}

void inline symb_put_simple_int_0(symb_coder_t *c, int val, int min, int max, uint64_t *count) {
  assert(min <= max);
  assert(val >= min);
  assert(val <= max);
  int x = val-min;
  if (x == 0) {
        symb_put_simple_bit(c, 0, count);
  } else {
        symb_put_simple_bit(c, 1, count);
        symb_put_simple_int(c, x, 1, max-min, count);
  }        
}

int inline symb_get_simple_int_0(symb_coder_t *c, int min, int max) {
  assert(min <= max);
  int x = 0;
  int nonzero = symb_get_simple_bit(c);
  if (nonzero) x = symb_get_simple_int(c,1,max-min);
  return x+min;
}



// do not output/input anything, but update chance tables anyway
void inline symb_skip_bit(symb_coder_t *c, chs_t *chs, int val) {
  chs_put(chs,&c->table,val);
}

void symb_cp(symb_chs_t *from,symb_chs_t *to) {
   for (int i=0; i<30; i++) {
      chs_cp(&from->chs[i],&to->chs[i]);
   }
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
  for (int i = 9; i < 30; i++) {
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
//static uint64_t pos = 0, neg = 0, neut = 0;
//static int64_t nb_diffs = 0;
//, tot_diffs = 0, tot_sq_diffs = 0;

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

/*
// alternative method: just write the number in binary
inline void symb_put_int_limited_v_a(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, int virtua, colors *cols, int channel, pixel_t *pixel, int qz, int guess) {
//  fprintf(stderr,"jif_put_symb: %4i in [%4i..%4i]\n",*val,min,max);
  nb_diffs++;
  int result=reduce_range(cols,&min,&max,channel,pixel,qz,guess);
  if (result == 0) fprintf(stderr,"PROBLEM output: trying to output a color that does not exist (channel %i, qz=%i)\n",channel,qz);
  if (result == 2) { if (*val == min) { return; } else { fprintf(stderr,"PROBLEM output: color is known but not correct!\n");  }   }
  assert(min<max);
  assert(*val>=min);
  assert(*val<=max);

  int tmin = 0;  int tmax = 0;
  int canbezero=reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  if (!(max >= 0 && min <= 0)) canbezero = 0;

  if (canbezero) {
      symb_put_bit_v(c,&sc->chs[0],*val == 0,count,virtua);
      if (*val == 0) return;
  }
  assert(*val != 0);
  tmin = 1;
  tmax = max;
  int canbepos=reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  tmin = min;
  tmax = -1;
  int canbeneg=reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  if (!canbepos && !canbeneg) return;   // has to be zero
  int negative = (*val < 0);
  if (canbepos && canbeneg && min*max<0) {
       symb_put_bit_v(c,&sc->chs[1],negative,count,virtua);
  }
  const int a = abs(*val);
  if (negative) if (max > -1) max = -1;
  if (!negative) if (min < 1) min = 1;
  int amax = abs(max) > abs(min) ? abs(max) : abs(min);
//  int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
  int bmax = log2_tab[amax];
//  int bmin = log2_tab[amin];

  int run = 0;
  int left = (1 << (bmax+1)) - 1;
  for (int i = bmax; i>=0; i--) {
      int minval = (negative ? -(run+left) : run);
      int maxval = (negative ? -run : run+left);
      left ^= (1 << i);
      if (min > minval) minval = min;
      if (max < maxval) maxval = max;
      int rr=reduce_range(cols,&minval,&maxval,channel,pixel,qz,guess);
      if (rr == 2) { run = (negative ? -minval : minval); break;}
      if (rr == 0) { fprintf(stderr,"PROBLEM output: impossible mantissa\n"); assert(1>2);}
      int bit = 1;
      int tmin1 = (negative ? minval : run + (1 << i));
      int tmax1 = (negative ? -(run + (1 << i)) : maxval);
      int tmin0 = (negative ? -(run+left) : minval);
      int tmax0 = (negative ? maxval : run + left);
      if (! reduce_range(cols,&tmin1,&tmax1,channel,pixel,qz,guess)) { // 1-bit would cause overflow
        bit = 0;
        symb_skip_bit(c,&sc->chs[2 + i],0);
      } else if (! reduce_range(cols,&tmin0,&tmax0,channel,pixel,qz,guess)) { // 0-bit would cause underflow
        symb_skip_bit(c,&sc->chs[2 + i],1);
      } else { // both 0 and 1 are possible
        if (i >= bmax-1-nb_bits) {
          bit = (a >> i) & 1;
          symb_put_bit_v(c,&sc->chs[2 + i],bit,count,virtua);
        } else if (i == bmax-nb_bits-1) {
          bit = 1;
        } else {
          bit = 0;
        }
      }
      run |= (bit << i);
  }
  if (*val != (negative?-run:run)) {
      fprintf(stderr,"WARNING: orig=%i, outputted=%i, range=%i..%i, a=%i, e=%i\n",*val,(*val<0?-run:run),min,max,a,bmax);
  }
  *val = clamp((negative?-run:run),min,max);

  assert(*val>=min);
  assert(*val<=max);
}




inline int symb_get_int_limited_c_a(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits, uint64_t *count, colors *cols, int channel, pixel_t *pixel, int qz, int guess) {
//  fprintf(stderr,"jif_get_symb: ???? in [%4i..%4i]\n",min,max);
  nb_diffs++;
  int result=reduce_range(cols,&min,&max,channel,pixel,qz,guess);
  if (result == 0) fprintf(stderr,"PROBLEM input: trying to output a color that does not exist (channel %i, qz=%i)\n",channel,qz);
  if (result == 2) { return min; }
  if (min == max) return min;
  assert(min<max);

  int tmin = 0;  int tmax = 0;
  int canbezero=reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  if (!(max >= 0 && min <= 0)) canbezero = 0;

  if (canbezero) {
      if (symb_get_bit_c(c,&sc->chs[0],count)) return 0;
  }

  tmin = 1;
  tmax = max;
  int canbepos=reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  tmin = min;
  tmax = -1;
  int canbeneg=reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  if (!canbepos && !canbeneg) return 0;
  int negative = 0;
  if (canbepos && canbeneg && min*max<0) {
       negative = (symb_get_bit_c(c,&sc->chs[1],count));
  } else {
       if (min>0 || !canbeneg) negative = 0;
       if (max<0 || !canbepos) negative = 1;
  }

  if (negative) if (max > -1) max = -1;
  if (!negative) if (min < 1) min = 1;
  int amax = abs(max) > abs(min) ? abs(max) : abs(min);
//  int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
  int bmax = log2_tab[amax];
//  int bmin = log2_tab[amin];

  int run = 0;
  int left = (1 << (bmax+1)) - 1;
  for (int i = bmax; i>=0; i--) {
      int minval = (negative ? -(run+left) : run);
      int maxval = (negative ? -run : run+left);
      left ^= (1 << i);
      if (min > minval) minval = min;
      if (max < maxval) maxval = max;
      int rr=reduce_range(cols,&minval,&maxval,channel,pixel,qz,guess);
      if (rr == 2) { run = (negative ? -minval : minval); break;}
      if (rr == 0) { fprintf(stderr,"PROBLEM input: impossible mantissa\n"); assert(1>2);}
      int bit = 1;
      int tmin1 = (negative ? minval : run + (1 << i));
      int tmax1 = (negative ? -(run + (1 << i)) : maxval);
      int tmin0 = (negative ? -(run+left) : minval);
      int tmax0 = (negative ? maxval : run + left);
      if (! reduce_range(cols,&tmin1,&tmax1,channel,pixel,qz,guess)) { // 1-bit would cause overflow
        bit = 0;
        symb_skip_bit(c,&sc->chs[2 + i],0);
      } else if (! reduce_range(cols,&tmin0,&tmax0,channel,pixel,qz,guess)) { // 0-bit would cause underflow
        symb_skip_bit(c,&sc->chs[2 + i],1);
      } else { // both 0 and 1 are possible
        if (i >= bmax-1-nb_bits) {
          bit = symb_get_bit_c(c,&sc->chs[2 + i],count);
        } else if (i == bmax-nb_bits-1) {
          bit = 1;
        } else {
          bit = 0;
        }
      }
      run |= (bit << i);
  }

    int ret = clamp((negative?-run:run),min,max);

    // output
    if (ret<min || ret>max) fprintf(stderr,"jif_get_symb: %4i in [%4i..%4i]\n",ret,min,max);
    assert(ret>=min);
    assert(ret<=max);
    return ret;

}

*/





/* write out a whole symbol to the range encoder using exponent/mantissa/sign representation
 - a zero bit is written if zero is within [min..max]
 - exponent bits are written in unary notation as long as no overflow is implied by them
 - mantissa bits are written only when necessary
 - a sign bit is written if both val and -val are within [min..max]
 */

static inline void symb_put_int_limited_v_a(symb_coder_t *c, symb_chs_t *sc, int *val, const int min, const int max, uint64_t *count, const int nb_bits, const int virtua, bitvector *bv) {
//colors *cols, int channel, pixel_t *pixel, int qz, int guess) {
//  fprintf(stderr,"jif_put_symb: %4i in [%4i..%4i]\n",*val,min,max);
  // initialization

  assert(min<=max);
  assert(*val>=min);
  assert(*val<=max);



//  int result=count_range(bv,min,max,min);
/*  if (result == 0) {
        fprintf(stderr,"PROBLEM output: trying to output a color that does not exist\n");
        return;
  }*/
/*  assert(result>0);
  if (result == 1) {
    assert(*val == get_first_in_range(bv,min,max,min));
    return;*/
/*
    if (*val == get_first_in_range(min,max,min)) {
//        fprintf(stderr,"YES: color range reduced to singleton\n");
//        virtua = 1;
           return;
    } else {
        fprintf(stderr,"PROBLEM output: color is known but not correct!\n");
        assert(0);
    }
*/
//  }


  if (*val) { // nonzero
    const unsigned int a = abs(*val);
    const unsigned int e = log2_tab[a];

    int canbezero=0;
    if ((max >= 0 && min <= 0)) {canbezero=check_position(bv,0);}

    // zeroness (chance state 0)
    if (canbezero) { // zero is possible
      symb_put_bit_v(c,&sc->chs[0],0,count,virtua);
//    } else {
//      symb_skip_bit(c,&sc->chs[0],0);
    }
//    if (!canbepos) if (max > -1) max = -1;
//    if (!canbeneg) if (min < 1) min = 1;
    unsigned int amax = abs(max) > abs(min) ? abs(max) : abs(min);
    unsigned int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
    unsigned int bmax = log2_tab[amax];
    unsigned int bmin = log2_tab[amin];


    // unary encoding of exponent (chance states 1..9)
    assert(e<=9);
    unsigned int mi = bmin;
    unsigned int again=1;

//    fprintf(stderr,"outputting exponent: bmin=%i, bmax=%i, e=%i\n",bmin,bmax,e);
    while (mi < bmax && again) {
      // check bigger exponents
      if (!check_range(bv,(1 << (mi+1)),(1 << (bmax+1))-1)   // bigger exponents, positive sign is imposible
        && !check_range(bv,-((1 << (bmax+1))-1),-(1 << (mi+1)))) {  // bigger exponents, negative sign is imposible
        assert(mi==e);
//        symb_skip_bit(c,&sc->chs[1 + mi],0);
        break;
      }
      // check current exponent
      if (!check_range(bv,(1 << (mi)),(1 << (mi))-1 + (1 << (mi)))          // this exponent is impossible with positive sign
         && !check_range(bv,-((1 << (mi))-1) -(1 << (mi)),-(1 << (mi)))) {  // this exponent is impossible with negative sign
        // this exponent is impossible, try next
 //       fprintf(stderr,"exponent cannot be %i\n",i);
//            symb_skip_bit(c,&sc->chs[1 + mi],1);
      } else {
            if (mi == e) again=0;
            symb_put_bit_v(c,&sc->chs[1+mi],again,count,virtua);
      }
      if (again) mi++;
    }
    assert(mi == e);
//    int rp=count_range(bv,(1 << e),(1 << e)-1 + (1 << e),min);
//    int rn=count_range(bv,-((1 << e)-1) -(1 << e),-(1 << e),min);
//    if (mi != e) fprintf(stderr,"PROBLEM output: exponent: bmin=%i, bmax=%i, expected=%i, real=%i\n",bmin,bmax,mi,e);
//    if (rp == 0 && rn == 0) fprintf(stderr,"PROBLEM output: impossible exponent: e=%i\n",e);

//    if (e < bmax) symb_put_bit_v(c,&sc->chs[1 + i],0,count,virtua);
//    fprintf(stderr,"exponent: %i\n",e);
    int run = (1 << e);
    int negative = (*val < 0);
 //    int sign = (negative ? -1 : 1);

//    rp = 1; rn = 1;
    if (min*max<0 && run<=max && -run>=min 
    && check_range(bv,(1 << e),(1 << e)-1 + (1 << e))   
    && check_range(bv,-((1 << e)-1) -(1 << e),-(1 << e)) ) {   
      // otherwise one option is impossible, no need to output sign
     // sign (chance states 11..20)
       symb_put_bit_v(c,&sc->chs[11 + e],negative,count,virtua);
    } //else { symb_skip_bit(c,&sc->chs[11+e],negative); }
    
//    if (negative && rn == 1) return;
//    if (!negative && rp == 1) return;
    
    // mantissa (chance states 21..29)
    int left = (1 << e) - 1;

    for (unsigned int i = e; i>0; ) {
      int bit = 1;
      int minval = (negative ? -(run+left) : run);
      int maxval = (negative ? -run : run+left);
      left ^= (1 << (--i));
//      fprintf(stderr,"jif_put_symb3: %4i in [%4i..%4i]\n",*val,minval,maxval);
      if (min > minval) minval = min;
      if (max < maxval) maxval = max;
      int rr=count_range(bv,minval,maxval);
//      fprintf(stderr,"jif_put_symb4: %4i in [%4i..%4i]\n",*val,minval,maxval);
      if (rr == 1) { assert(*val == get_first_in_range(bv,minval,maxval)); return; } //*val = get_first_in_range(bv,minval,maxval,min); return; }
//      if (rr == 0) {fprintf(stderr,"PROBLEM output: impossible mantissa\n"); assert(1>2);}
      assert(rr > 1);
      int tmin1 = (negative ? minval : run + (1 << i));
      int tmax1 = (negative ? -(run + (1 << i)) : maxval);
      int tmin0 = (negative ? -(run+left) : minval);
      int tmax0 = (negative ? maxval : run + left);
      if (! check_range(bv,tmin1,tmax1)) { // 1-bit would cause overflow
        bit = 0;
//        symb_skip_bit(c,&sc->chs[21 + i],0);
      } else if (! check_range(bv,tmin0,tmax0)) { // 0-bit would cause underflow
//        symb_skip_bit(c,&sc->chs[21 + i],1);
      } else { // both 0 and 1 are possible
//        if (i+nb_bits >= e) {
          bit = (a >> i) & 1;
          symb_put_bit_v(c,&sc->chs[21 + i],bit,count,virtua);
/*        } else if (i+nb_bits+1 == e) {
          bit = 1;
        } else {
          bit = 0;
        }*/
      }
      run |= (bit << i);
    }

/*    if (*val != (*val<0?-run:run)) {
        fprintf(stderr,"WARNING: orig=%i, outputted=%i, range=%i..%i, a=%i, e=%i\n",*val,(*val<0?-run:run),min,max,a,e);
//        assert(5<3);
    }
*/
    assert(    *val == (*val<0?-run:run) );
//    *val = clamp((*val<0?-run:run),min,max);


  } else { // zero
    if (check_range(bv,min,-1) || check_range(bv,1,max) )       // non-zero is possible
            symb_put_bit_v(c,&sc->chs[0],1,count,virtua);
//    if (!canbezero) fprintf(stderr,"PROBLEM output: was zero and cannot be\n");
//    assert(canbezero);
  }
  assert(*val>=min);
  assert(*val<=max);
}

// read in a whole symbol from the range encoder using exponent/mantissa/sign representation
inline int symb_get_int_limited_c_a(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits, uint64_t *count, bitvector *bv) {
//colors *cols, int channel, pixel_t *pixel, int qz, int guess) {

//  fprintf(stderr,"jif_get_symb: ???? in [%4i..%4i]\n",min,max);
//  assert(min<=max);


//  int result=count_range(bv,min,max,min);
/*  if (result == 0) {
        fprintf(stderr,"PROBLEM input: trying to input a color that does not exist (channel %i)\n",channel);
        return (min+max)/2;
  }*/
/*  assert(result>0);
  if (result == 1) {
    return get_first_in_range(bv,min,max,min);
  }*/
//  fprintf(stderr,"jif_get_symb2: ???? in [%4i..%4i]\n",min,max);
  // initialization
  assert(min<=max);
  if (min == max) return min;

  assert(check_range(bv,min,max));

  int canbezero=0;
  if ((max >= 0 && min <= 0)) {canbezero=check_position(bv,0);}

/*
  tmin = 1;
  tmax = max;
  int canbepos=1;//reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  tmin = min;
  tmax = -1;
  int canbeneg=1;//reduce_range(cols,&tmin,&tmax,channel,pixel,qz,guess);
  if (!canbepos && !canbeneg) return 0;   // has to be zero
*/



  if (canbezero && 
        ( (!check_range(bv,min,-1) && !check_range(bv,1,max))   // must be zero
        || symb_get_bit_c(c,&sc->chs[0],count)) ) {                       // is zero
//    fprintf(stderr,"jif_get_symb:    0 in [%4i..%4i]\n",min,max);
    return 0;
  } else { // nonzero
//      if (!canbezero) symb_skip_bit(c,&sc->chs[0],0);

    // unary encoding of exponent (chance states 1..10)

//    if (!canbepos) if (max > -1) max = -1;
//    if (!canbeneg) if (min < 1) min = 1;
    unsigned int amax = abs(max) > abs(min) ? abs(max) : abs(min);
    unsigned int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
    unsigned int bmax = log2_tab[amax];
    unsigned int bmin = log2_tab[amin];

    unsigned int e = bmin;
    unsigned int again=1;
    while (e < bmax && again) {
      // check bigger exponents
      if (!check_range(bv,(1 << (e+1)),(1 << (bmax+1))-1)   // bigger exponents, positive sign is imposible
        && !check_range(bv,-((1 << (bmax+1))-1),-(1 << (e+1)))) {  // bigger exponents, negative sign is imposible
//        symb_skip_bit(c,&sc->chs[1 + e],0);
        break;
      }
      // check current exponent
      if (!check_range(bv,(1 << (e)),(1 << (e))-1 + (1 << (e)))          // this exponent is impossible with positive sign
         && !check_range(bv,-((1 << (e))-1) -(1 << (e)),-(1 << (e)))) {  // this exponent is impossible with negative sign
        // this exponent is impossible, try next
//        symb_skip_bit(c,&sc->chs[1 + e],1);
      } else {
        again=symb_get_bit_c(c,&sc->chs[1+e],count);
      }
      if (again) e++;
    }

    unsigned int rp = 1;
    unsigned int rn = 1;



    assert(e<=9);
//    fprintf(stderr,"exponent: %i\n",e);
    //unsigned 
    int run = (1 << e);
//    rp = 1; rn = 1;
    int sign = 0;
    if ( min*max < 0 && run<=max && -run>=min
    && (rp=check_range(bv,(1 << e),(1 << e)-1 + (1 << e))) > 0 
    && (rn=check_range(bv,-((1 << e)-1) -(1 << e),-(1 << e))) > 0) {     // otherwise one option is impossible, no need to input sign
      // sign (chance states 11..20)
        if (symb_get_bit_c(c,&sc->chs[11 + e],count)) {
          sign = -1;
        } else {
          sign = 1;
        }
    } else {
        sign = 1;
        if (run>max) sign = -1;
        if (rp == 0) sign = -1;
        if (rn == 0) sign = 1;
        assert(!(rp==0 && rn ==0));
        //if (rp == 0 && rn == 0) fprintf(stderr,"PROBLEM input: impossible sign\n");
//        symb_skip_bit(c,&sc->chs[11+e],(sign < 0));
    }


    int negative = (sign < 0);

//    if (negative && rn == 1) return get_first_in_range(bv,emin2,emax2,min);
//    if (!negative && rp == 1) return get_first_in_range(bv,emin,emax,min);


    // mantissa (chance states 21..29)
    int left = (1 << e) - 1;
    for (unsigned int i = e; i>0; ) {
      int minval = (negative ? -(run+left) : run);
      int maxval = (negative ? -run : run+left);
      left ^= (1 << (--i));
      if (min > minval) minval = min;
      if (max < maxval) maxval = max;
      unsigned int rr=count_range(bv,minval,maxval);
      if (rr == 1) { return get_first_in_range(bv,minval,maxval); }
//      if (rr == 0) { fprintf(stderr,"PROBLEM input: impossible mantissa\n"); assert(1>2);}
      assert(rr>1);
      int tmin1 = (negative ? minval : run + (1 << i));
      int tmax1 = (negative ? -(run + (1 << i)) : maxval);
      int tmin0 = (negative ? -(run+left) : minval);
      int tmax0 = (negative ? maxval : run + left);
      int bit = 1;
      if (! check_range(bv,tmin1,tmax1)) { // 1-bit would cause overflow
        bit = 0;
//        symb_skip_bit(c,&sc->chs[21 + i],0);
//      } else if (run + left < amin) { // 0-bit would cause underflow
      } else if (! check_range(bv,tmin0,tmax0)) { // 0-bit would cause underflow
//        symb_skip_bit(c,&sc->chs[21 + i],1);
      } else { // both 0 and 1 are possible
//        if (i+nb_bits >= e) {
          bit = symb_get_bit_c(c,&sc->chs[21 + i],count);
/*        } else if (i+nb_bits+1 == e) {
          bit = 1;
        } else {
          bit = 0;
        }  */
      }
      run |= (bit << i);
    }
    int ret = clamp(run*sign,min,max);

    // output
//    if (ret<min || ret>max) fprintf(stderr,"jif_get_symb: %4i in [%4i..%4i]\n",ret,min,max);
    assert(ret>=min);
    assert(ret<=max);
    return ret;
  }
}





inline void symb_put_int_limited_v(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, int virtua) {
//   disable_vector();
   symb_put_int_limited_v_a(c, sc, val, min, max, count, nb_bits, virtua, NULL);
//   enable_vector();
   //NULL, -1, NULL, 1, 0);
}
inline void symb_put_int_limited_a(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, bitvector *bv) {
//   symb_put_int_limited_v_a(c, sc, val, min, max, count, nb_bits, 0, bv);

  if (*val) { // nonzero
    const unsigned int a = abs(*val);
    const unsigned int e = log2_tab[a];

    int canbezero=0;
    if ((max >= 0 && min <= 0)) {canbezero=check_position(bv,0);}

    // zeroness (chance state 0)
    if (canbezero) { // zero is possible
      symb_put_bit_v(c,&sc->chs[0],0,count,0);
//    } else {
//      symb_skip_bit(c,&sc->chs[0],0);
    }
//    if (!canbepos) if (max > -1) max = -1;
//    if (!canbeneg) if (min < 1) min = 1;
    unsigned int amax = abs(max) > abs(min) ? abs(max) : abs(min);
    unsigned int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
    unsigned int bmax = log2_tab[amax];
    unsigned int bmin = log2_tab[amin];


    // unary encoding of exponent (chance states 1..9)
    assert(e<=9);
    unsigned int mi = bmin;
    unsigned int again=1;

//    fprintf(stderr,"outputting exponent: bmin=%i, bmax=%i, e=%i\n",bmin,bmax,e);
    while (mi < bmax && again) {
      // check bigger exponents
      if (!check_range(bv,(1 << (mi+1)),(1 << (bmax+1))-1)   // bigger exponents, positive sign is imposible
        && !check_range(bv,-((1 << (bmax+1))-1),-(1 << (mi+1)))) {  // bigger exponents, negative sign is imposible
        assert(mi==e);
//        symb_skip_bit(c,&sc->chs[1 + mi],0);
        break;
      }
      // check current exponent
      if (!check_range(bv,(1 << (mi)),(1 << (mi))-1 + (1 << (mi)))          // this exponent is impossible with positive sign
         && !check_range(bv,-((1 << (mi))-1) -(1 << (mi)),-(1 << (mi)))) {  // this exponent is impossible with negative sign
        // this exponent is impossible, try next
 //       fprintf(stderr,"exponent cannot be %i\n",i);
//            symb_skip_bit(c,&sc->chs[1 + mi],1);
      } else {
            if (mi == e) again=0;
            symb_put_bit_v(c,&sc->chs[1+mi],again,count,0);
      }
      if (again) mi++;
    }
    assert(mi == e);
//    int rp=count_range(bv,(1 << e),(1 << e)-1 + (1 << e),min);
//    int rn=count_range(bv,-((1 << e)-1) -(1 << e),-(1 << e),min);
//    if (mi != e) fprintf(stderr,"PROBLEM output: exponent: bmin=%i, bmax=%i, expected=%i, real=%i\n",bmin,bmax,mi,e);
//    if (rp == 0 && rn == 0) fprintf(stderr,"PROBLEM output: impossible exponent: e=%i\n",e);

//    if (e < bmax) symb_put_bit_v(c,&sc->chs[1 + i],0,count,0);
//    fprintf(stderr,"exponent: %i\n",e);
    int run = (1 << e);
    int negative = (*val < 0);
 //    int sign = (negative ? -1 : 1);

//    rp = 1; rn = 1;
    if (min*max<0 && run<=max && -run>=min 
    && check_range(bv,(1 << e),(1 << e)-1 + (1 << e))   
    && check_range(bv,-((1 << e)-1) -(1 << e),-(1 << e)) ) {   
      // otherwise one option is impossible, no need to output sign
     // sign (chance states 11..20)
       symb_put_bit_v(c,&sc->chs[11 + e],negative,count,0);
    } //else { symb_skip_bit(c,&sc->chs[11+e],negative); }
    
//    if (negative && rn == 1) return;
//    if (!negative && rp == 1) return;
    
    // mantissa (chance states 21..29)
    int left = (1 << e) - 1;

    for (unsigned int i = e; i>0; ) {
      int bit = 1;
      int minval = (negative ? -(run+left) : run);
      int maxval = (negative ? -run : run+left);
      left ^= (1 << (--i));
//      fprintf(stderr,"jif_put_symb3: %4i in [%4i..%4i]\n",*val,minval,maxval);
      if (min > minval) minval = min;
      if (max < maxval) maxval = max;
      int rr=count_range(bv,minval,maxval);
//      fprintf(stderr,"jif_put_symb4: %4i in [%4i..%4i]\n",*val,minval,maxval);
      if (rr == 1) { assert(*val == get_first_in_range(bv,minval,maxval)); return; } //*val = get_first_in_range(bv,minval,maxval,min); return; }
//      if (rr == 0) {fprintf(stderr,"PROBLEM output: impossible mantissa\n"); assert(1>2);}
      assert(rr > 1);
      int tmin1 = (negative ? minval : run + (1 << i));
      int tmax1 = (negative ? -(run + (1 << i)) : maxval);
      int tmin0 = (negative ? -(run+left) : minval);
      int tmax0 = (negative ? maxval : run + left);
      if (! check_range(bv,tmin1,tmax1)) { // 1-bit would cause overflow
        bit = 0;
//        symb_skip_bit(c,&sc->chs[21 + i],0);
      } else if (! check_range(bv,tmin0,tmax0)) { // 0-bit would cause underflow
//        symb_skip_bit(c,&sc->chs[21 + i],1);
      } else { // both 0 and 1 are possible
//        if (i+nb_bits >= e) {
          bit = (a >> i) & 1;
          symb_put_bit_v(c,&sc->chs[21 + i],bit,count,0);
/*        } else if (i+nb_bits+1 == e) {
          bit = 1;
        } else {
          bit = 0;
        }*/
      }
      run |= (bit << i);
    }

/*    if (*val != (*val<0?-run:run)) {
        fprintf(stderr,"WARNING: orig=%i, outputted=%i, range=%i..%i, a=%i, e=%i\n",*val,(*val<0?-run:run),min,max,a,e);
//        assert(5<3);
    }
*/
    assert(    *val == (*val<0?-run:run) );
//    *val = clamp((*val<0?-run:run),min,max);


  } else { // zero
    if (check_range(bv,min,-1) || check_range(bv,1,max) )       // non-zero is possible
            symb_put_bit_v(c,&sc->chs[0],1,count,0);
//    if (!canbezero) fprintf(stderr,"PROBLEM output: was zero and cannot be\n");
//    assert(canbezero);
  }
  assert(*val>=min);
  assert(*val<=max);

}
inline void symb_put_int_limited_a1(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits, bitvector *bv) {
//   symb_put_int_limited_v_a(c, sc, val, min, max, count, nb_bits, 1, bv);

  if (*val) { // nonzero
    const unsigned int a = abs(*val);
    const unsigned int e = log2_tab[a];

    int canbezero=0;
    if ((max >= 0 && min <= 0)) {canbezero=check_position(bv,0);}

    // zeroness (chance state 0)
    if (canbezero) { // zero is possible
      symb_put_bit_v(c,&sc->chs[0],0,count,1);
//    } else {
//      symb_skip_bit(c,&sc->chs[0],0);
    }
//    if (!canbepos) if (max > -1) max = -1;
//    if (!canbeneg) if (min < 1) min = 1;
    unsigned int amax = abs(max) > abs(min) ? abs(max) : abs(min);
    unsigned int amin = (max >= 0 && min <= 0) ? 0 : (abs(max) < abs(min) ? abs(max) : abs(min));
    unsigned int bmax = log2_tab[amax];
    unsigned int bmin = log2_tab[amin];


    // unary encoding of exponent (chance states 1..9)
    assert(e<=9);
    unsigned int mi = bmin;
    unsigned int again=1;

//    fprintf(stderr,"outputting exponent: bmin=%i, bmax=%i, e=%i\n",bmin,bmax,e);
    while (mi < bmax && again) {
      // check bigger exponents
      if (!check_range(bv,(1 << (mi+1)),(1 << (bmax+1))-1)   // bigger exponents, positive sign is imposible
        && !check_range(bv,-((1 << (bmax+1))-1),-(1 << (mi+1)))) {  // bigger exponents, negative sign is imposible
        assert(mi==e);
//        symb_skip_bit(c,&sc->chs[1 + mi],0);
        break;
      }
      // check current exponent
      if (!check_range(bv,(1 << (mi)),(1 << (mi))-1 + (1 << (mi)))          // this exponent is impossible with positive sign
         && !check_range(bv,-((1 << (mi))-1) -(1 << (mi)),-(1 << (mi)))) {  // this exponent is impossible with negative sign
        // this exponent is impossible, try next
 //       fprintf(stderr,"exponent cannot be %i\n",i);
//            symb_skip_bit(c,&sc->chs[1 + mi],1);
      } else {
            if (mi == e) again=0;
            symb_put_bit_v(c,&sc->chs[1+mi],again,count,1);
      }
      if (again) mi++;
    }
    assert(mi == e);
//    int rp=count_range(bv,(1 << e),(1 << e)-1 + (1 << e),min);
//    int rn=count_range(bv,-((1 << e)-1) -(1 << e),-(1 << e),min);
//    if (mi != e) fprintf(stderr,"PROBLEM output: exponent: bmin=%i, bmax=%i, expected=%i, real=%i\n",bmin,bmax,mi,e);
//    if (rp == 0 && rn == 0) fprintf(stderr,"PROBLEM output: impossible exponent: e=%i\n",e);

//    if (e < bmax) symb_put_bit_v(c,&sc->chs[1 + i],0,count,virtua);
//    fprintf(stderr,"exponent: %i\n",e);
    int run = (1 << e);
    int negative = (*val < 0);
 //    int sign = (negative ? -1 : 1);

//    rp = 1; rn = 1;
    if (min*max<0 && run<=max && -run>=min 
    && check_range(bv,(1 << e),(1 << e)-1 + (1 << e))   
    && check_range(bv,-((1 << e)-1) -(1 << e),-(1 << e)) ) {   
      // otherwise one option is impossible, no need to output sign
     // sign (chance states 11..20)
       symb_put_bit_v(c,&sc->chs[11 + e],negative,count,1);
    } //else { symb_skip_bit(c,&sc->chs[11+e],negative); }
    
//    if (negative && rn == 1) return;
//    if (!negative && rp == 1) return;
    
    // mantissa (chance states 21..29)
    int left = (1 << e) - 1;

    for (unsigned int i = e; i>0; ) {
      int bit = 1;
      int minval = (negative ? -(run+left) : run);
      int maxval = (negative ? -run : run+left);
      left ^= (1 << (--i));
//      fprintf(stderr,"jif_put_symb3: %4i in [%4i..%4i]\n",*val,minval,maxval);
      if (min > minval) minval = min;
      if (max < maxval) maxval = max;
      int rr=count_range(bv,minval,maxval);
//      fprintf(stderr,"jif_put_symb4: %4i in [%4i..%4i]\n",*val,minval,maxval);
      if (rr == 1) { assert(*val == get_first_in_range(bv,minval,maxval)); return; } //*val = get_first_in_range(bv,minval,maxval,min); return; }
//      if (rr == 0) {fprintf(stderr,"PROBLEM output: impossible mantissa\n"); assert(1>2);}
      assert(rr > 1);
      int tmin1 = (negative ? minval : run + (1 << i));
      int tmax1 = (negative ? -(run + (1 << i)) : maxval);
      int tmin0 = (negative ? -(run+left) : minval);
      int tmax0 = (negative ? maxval : run + left);
      if (! check_range(bv,tmin1,tmax1)) { // 1-bit would cause overflow
        bit = 0;
//        symb_skip_bit(c,&sc->chs[21 + i],0);
      } else if (! check_range(bv,tmin0,tmax0)) { // 0-bit would cause underflow
//        symb_skip_bit(c,&sc->chs[21 + i],1);
      } else { // both 0 and 1 are possible
//        if (i+nb_bits >= e) {
          bit = (a >> i) & 1;
          symb_put_bit_v(c,&sc->chs[21 + i],bit,count,1);
/*        } else if (i+nb_bits+1 == e) {
          bit = 1;
        } else {
          bit = 0;
        }*/
      }
      run |= (bit << i);
    }

/*    if (*val != (*val<0?-run:run)) {
        fprintf(stderr,"WARNING: orig=%i, outputted=%i, range=%i..%i, a=%i, e=%i\n",*val,(*val<0?-run:run),min,max,a,e);
//        assert(5<3);
    }
*/
    assert(    *val == (*val<0?-run:run) );
//    *val = clamp((*val<0?-run:run),min,max);


  } else { // zero
    if (check_range(bv,min,-1) || check_range(bv,1,max) )       // non-zero is possible
            symb_put_bit_v(c,&sc->chs[0],1,count,1);
//    if (!canbezero) fprintf(stderr,"PROBLEM output: was zero and cannot be\n");
//    assert(canbezero);
  }
  assert(*val>=min);
  assert(*val<=max);

}
inline void symb_put_int_limited(symb_coder_t *c, symb_chs_t *sc, int *val, int min, int max, uint64_t *count, int nb_bits) {
   symb_put_int_limited_v(c, sc, val, min, max, count, nb_bits, 0);
}
void symb_put_int(symb_coder_t *c, symb_chs_t *sc, int val, int min, int max, uint64_t *count) {
//   fprintf(stderr,"putting %i\n",val);
      assert(max<1023);
      symb_put_int_limited(c,sc, &val, min, max, count, 64);
//   fprintf(stderr,"Sput %i\n",val);
}

void symb_put_big_int(symb_coder_t *c, symb_chs_t *sc1, symb_chs_t *sc2, int val, int min, int max, uint64_t *count) {
      int rest = val % 1000;
      symb_put_int_limited(c, sc1, &rest, 0, 999, count, 64);
      int div = val / 1000;
      symb_put_int_limited(c, sc2, &div, min/1000, max/1000, count, 64);
      assert ( val == 1000*div+rest );
//   fprintf(stderr,"put %i = 1000*%i + %i\n",val,div,rest);
}


inline int symb_get_int_limited_c(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits, uint64_t *count) {
// disable_vector();
 return symb_get_int_limited_c_a(c, sc, min, max, nb_bits, count, NULL);
// enable_vector();
 //NULL, -1, NULL, 1, 0);
}
inline int symb_get_int_limited(symb_coder_t *c, symb_chs_t *sc, int min, int max, int nb_bits) {
 return symb_get_int_limited_c(c, sc, min, max, nb_bits, NULL);
}
int symb_get_int(symb_coder_t *c, symb_chs_t *sc, int min, int max) {
 int val=symb_get_int_limited(c, sc, min, max, 64);
// fprintf(stderr,"Sgot %i\n",val);
 return val;

}
int symb_get_big_int(symb_coder_t *c, symb_chs_t *sc1, symb_chs_t *sc2, int min, int max) {

// fprintf(stderr,"trying to get number between %i..%i\n",min,max);
    int rest = symb_get_int_limited(c, sc1, 0, 999, 64);
    int div = symb_get_int_limited(c, sc2, min/1000, max/1000, 64);
    int val = 1000*div+rest;
// fprintf(stderr,"got %i = 1000*%i + %i\n",val,div,rest);
 return val;

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
  fprintf(f,"Total bytes read/written (after header): %i\n",nbytes);
  chs_show_stats(f,&c->table);
}
