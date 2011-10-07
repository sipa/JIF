#ifndef _RAC_SYMBOL_H_
#define _RAC_SYMBOL_H_ 1

#include <stdio.h>

#include <vector>
#include <assert.h>
#include "util.h"
#include "chance.h"

typedef enum {
  BIT_ZERO,
  BIT_SIGN,
  BIT_EXP,
  BIT_MANT
} SymbolChanceBitType;

template <typename BitChance> class SymbolChance {
  int bits;
  std::vector<BitChance> chances;

public:

  BitChance inline &bitZero()      { return chances[0]; }
  BitChance inline &bitSign()      { return chances[1]; }

  // all exp bits 0         -> int(log2(val)) == 0  [ val == 1 ]
  // exp bits up to i are 1 -> int(log2(val)) == i+1
  BitChance inline &bitExp(int i)  { assert(i >= 0 && i < bits-1); return chances[2+i]; }
  BitChance inline &bitMant(int i) { assert(i >= 0 && i < bits); return chances[bits+1+i]; }

  BitChance inline &bit(SymbolChanceBitType typ, int i = 0) {
    switch (typ) {
      case BIT_ZERO: return bitZero();
      case BIT_SIGN: return bitSign();
      case BIT_EXP:  return bitExp(i);
      case BIT_MANT: return bitMant(i);
    }
  }

  SymbolChance(int bitsin) : bits(bitsin), chances(2*bits+1) { }
};

template <typename BitChance, typename RAC> class SimpleSymbolCoder {
private:
  SymbolChance<BitChance> ctx;
  RAC *rac;

public:
  SimpleSymbolCoder(RAC& racIn, int nBits) : ctx(nBits) {
    rac = &racIn;
  }

  void inline write(bool bit, SymbolChanceBitType typ, int i = 0) { 
    BitChance& bch = ctx.bit(typ,i);
    rac->write(bch.get(), bit);
    bch.put(bit);
  }

  bool inline read(SymbolChanceBitType typ, int i = 0) { 
    BitChance& bch = ctx.bit(typ,i);
    bool bit = rac->read(bch.get());
    bch.put(bit);
    return bit;
  }

};

int default_range_test(int a, int b) {
  if (b < a) return 0;
  assert(b >= a);
  return (b-a+1);
}

static inline int signed_test(int(*range_test)(int,int), int low, int high, int sign) {
  if (low > high) return 0;
  if (sign > 0) return range_test(low, high);
  return range_test(-high, -low);
}

template <typename SymbolCoder> int read_int(SymbolCoder& coder, int min, int max, int(*range_test)(int, int)) {
  assert(min<=max);
  assert(range_test(min,max));

  if (min == max) return min; // nodig?

  if (max >= 0 && min <= 0 && range_test(0,0)) {
    if (!signed_test(range_test,min,-1,1) && !signed_test(range_test,1,max,1)) return 0;
    if (coder.read(BIT_ZERO)) return 0;
  }

  bool sign = true;
  if (max > 0 && min < 0 && signed_test(range_test,min,-1,1) && signed_test(range_test,1,max,1)) {
    sign = coder.read(BIT_SIGN);
  } else {
    if (min < 0 || !signed_test(range_test,1,max,1)) sign = false;
  }
  if (sign && min <= 0) min = 1;
  if (!sign && max >= 0) max = -1;

  int emin = ilog2((sign ? abs(min) : abs(max)));
  int emax = ilog2((sign ? abs(max) : abs(min)));
//  fprintf(stdout,"min=%i,max=%i,emin:%i,emax:%i",min,max,emin,emax);
  int i = emin;
  while (i < emax) {
    // if exponent >i is impossible, we are done
    if (sign && !signed_test(range_test,1<<(i+1),max,1)) break;
    if (!sign && !signed_test(range_test,min,-(1<<(i+1)),1)) break;
    // if exponent i is possible, output the exponent bit
    if ((sign && range_test(1<<i,(1<<(i+1))-1)) || (!sign && range_test(1-(1<<(i+1)),-(1<<i)))) {
      if (coder.read(BIT_EXP,i)) break;
    }
    i++;
  }
  int e = i;
//  fprintf(stdout,"exp=%i",e);
  int have = (1 << e);
  int left = (1 << e)-1;
  for (int pos = e; pos>0;) {
    int bit = 1;
    int minval = (sign ? have : -(have+left));
    int maxval = (sign ? have+left : -have);
    if (min > minval) minval = min;
    if (max < maxval) maxval = max;
    left ^= (1 << (--pos));
    if (minval == maxval) return minval;
    int minval1 = (sign ? have+(1<<pos) : minval);
    int maxval1 = (sign ? maxval : -(have+(1<<pos)));
    int minval0 = (sign ? minval : -(have+left));
    int maxval0 = (sign ? have+left : maxval);
    if (!signed_test(range_test,minval1,maxval1,1)) { // 1-bit is impossible
      bit = 0;
    } else if (signed_test(range_test,minval0,maxval0,1)) { // 0-bit and 1-bit are both possible
      bit = coder.read(BIT_MANT,pos);
    }
    have |= (bit << pos);
  }
  int a = have;
  return (sign ? a : -a);
}

template <typename SymbolCoder> void write_int(SymbolCoder& coder, int min, int max, int(*range_test)(int, int), int &value) {
    assert(min<=max);
    assert(value>=min);
    assert(value<=max);
    assert(signed_test(range_test,min,max,1));

    // avoid doing anything if the value is already known (this line is optional; behavior should be identical if commented out)
    if (signed_test(range_test,min,max,1) == 1) return;

    if (value) { // value is nonzero
      // only output zero bit if value could also have been zero
      if (max >= 0 && min <= 0 && range_test(0,0)) coder.write(false,BIT_ZERO);
      int sign = (value > 0 ? 1 : 0);
      if (max > 0 && min < 0) {
        // only output sign bit if value can be both pos and neg
        if (range_test(min,-1) && range_test(1,max)) coder.write(sign,BIT_SIGN);
      }
      if (sign && min <= 0) min = 1;
      if (!sign && max >= 0) max = -1;
      const int a = abs(value);
      const int e = ilog2(a);
      int emin = ilog2((sign ? abs(min) : abs(max)));
      int emax = ilog2((sign ? abs(max) : abs(min)));
//      fprintf(stdout,"min=%i,max=%i,emin:%i,emax:%i",min,max,emin,emax);
      int i = emin;
      while (i < emax) {
        // if exponent >i is impossible, we are done
        if (sign && !range_test(1<<(i+1),max)) break;
        if (!sign && !range_test(min,-(1<<(i+1)))) break;
        // if exponent i is possible, output the exponent bit
        if ((sign && range_test(1<<i,(1<<(i+1))-1)) || (!sign && range_test(1-(1<<(i+1)),-(1<<i))))
            coder.write(i==e, BIT_EXP, i);
        if (i==e) break;
        i++;
      }
//      fprintf(stdout,"exp=%i",e);
      int have = (1 << e);
      int left = (1 << e)-1;
      for (int pos = e; pos>0;) {
        int bit = 1;
        int minabs = have, maxabs = have+left;
        left ^= (1 << (--pos));
        if (signed_test(range_test,minabs,maxabs,sign)==1) return;
        int minabs1 = have + (1<<pos), maxabs1 = maxabs;
        int minabs0 = minabs, maxabs0 = have + left;
        if (!signed_test(range_test,minabs1,maxabs1,sign)) { // 1-bit is impossible
           bit = 0;
        } else if (signed_test(range_test,minabs0,maxabs0,sign)) { // 0-bit and 1-bit are both possible
           bit = (a >> pos) & 1;
           coder.write(bit, BIT_MANT, pos);
        }
        have |= (bit << pos);
      }

    } else { // value is zero
      // only output zero bit if value could also have been nonzero
      if (range_test(min,-1) || range_test(1,max)) coder.write(true, BIT_ZERO);
    }
}


#endif
