#ifndef _RAC_SYMBOL_H_
#define _RAC_SYMBOL_H_ 1

#include <vector>
#include <assert.h>
#include "util.h"
#include "chance.h"

template <typename BitChance> class SymbolContext {
  int bits;
  std::vector<BitChance> chances;

public:
  BitChance inline &bitZero()      { return chances[0]; }
  BitChance inline &bitSign()      { return chances[1]; }

  // all exp bits 0         -> int(log2(val)) == 0  [ val == 1 ]
  // exp bits up to i are 1 -> int(log2(val)) == i+1
  BitChance inline &bitExp(int i)  { assert(i >= 0 && i < bits-1); return chances[2+i]; }
  BitChance inline &bitMant(int i) { assert(i >= 0 && i < bits); return chances[bits+1+i]; }

  SymbolContext(int bitsin) : bits(bitsin), chances(2*bits+1) { }
};

template <typename BitChance, typename RAC> class SimpleSymbolCoder {
private:
  SymbolContext<BitChance> ctx;
  RAC *rac;

  void inline write(BitChance& bch, bool bit) { rac->write(bch.get(), bit); bch.put(bit); }
  bool inline read(BitChance& bch) { bool bit = rac->read(bch.get()); bch.put(bit); return bit; }

public:
  SimpleSymbolCoder(RAC& racIn, int nBits) : ctx(nBits) {
    rac = &racIn;
  }

  void inline writeZero(bool bit) { write(ctx.bitZero(), bit); }
  void inline writeSign(bool bit) { write(ctx.bitSign(), bit); }
  void inline writeExp(int i, bool bit) { write(ctx.bitExp(i), bit); }
  void inline writeMant(int i, bool bit) { write(ctx.bitMant(i), bit); }
  bool inline readZero() { return read(ctx.bitZero()); }
  bool inline readSign() { return read(ctx.bitSign()); }
  bool inline readExp(int i) { return read(ctx.bitExp(i)); }
  bool inline readMant(int i) { return read(ctx.bitMant(i)); }
};

int default_range_test(int a, int b) {
  assert(b >= a);
  return (b-a+1);
}

template <typename SymbolCoder> int read_int(SymbolCoder& coder, int min, int max, int(*range_test)(int, int)) {
  assert(min<=max);
  assert(range_test(min,max));

  if (min == max) return min; // nodig?

  if (max >= 0 && min <= 0 && range_test(0,0)) {
    if (!range_test(min,-1) && !range_test(1,max)) return 0;
    if (coder.readZero()) return 0;
  }

  bool sign = true;
  if (max > 0 && min < 0 && range_test(min,-1) && range_test(1,max)) {
    sign = coder.readSign();
  } else {
    if (min < 0 || !range_test(1,max)) sign = false;
  }
  if (sign && min <= 0) min = 1;
  if (!sign && max >= 0) max = -1;
  
  unsigned int emin = ilog2((sign ? abs(min) : abs(max)));
  unsigned int emax = ilog2((sign ? abs(max) : abs(min)));
  unsigned int i = emin;
  while (i < emax) {
    // if exponent >i is impossible, we are done
    if (sign && !range_test(1<<(i+1),max)) break;
    if (!sign && !range_test(min,-(1<<(i+1)))) break;
    // if exponent i is possible, output the exponent bit
    if (sign && range_test(1<<i,1<<(i+1)-1) || !sign && !range_test(1-(1<<(i+1)),-(1<<i))) {
      if (coder.readExp(i)) break;
    }
    i++;
  }
  int e = i;
  int have = (1 << e);
  int left = (1 << e)-1;
  for (unsigned int pos = e; pos>0;) {
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
    if (!range_test(minval1,maxval1)) { // 1-bit is impossible
      bit = 0;
    } else if (range_test(minval0,maxval0)) { // 0-bit and 1-bit are both possible
      bit = coder.readMant(pos);
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
    assert(range_test(min,max));

    // avoid doing anything if the value is already known (this line is optional; behavior should be identical if commented out)
    if (range_test(min,max) == 1) return;

    if (value) { // value is nonzero
      // only output zero bit if value could also have been zero
      if (max >= 0 && min <= 0 && range_test(0,0)) coder.writeZero(false);
      bool sign = (value > 0 ? true : false);
      // only output sign bit if value can be both pos and neg
      if (max > 0 && min < 0 && range_test(min,-1) && range_test(1,max)) coder.writeSign(sign);
      if (sign && min <= 0) min = 1;
      if (!sign && max >= 0) max = -1;
      const unsigned int a = abs(value);
      const unsigned int e = ilog2(a);
      unsigned int emin = ilog2((sign ? abs(min) : abs(max)));
      unsigned int emax = ilog2((sign ? abs(max) : abs(min)));
      unsigned int i = emin;
      while (i < emax) {
        // if exponent >i is impossible, we are done
        if (sign && !range_test(1<<(i+1),max)) break;
        if (!sign && !range_test(min,-(1<<(i+1)))) break;
        // if exponent i is possible, output the exponent bit
        if ( (sign && range_test(1<<i,(1<<(i+1))-1))
         || (!sign && !range_test(1-(1<<(i+1)),-(1<<i))))
                coder.writeExp(i, i==e);
        if (i == e) break;
        i++;
      }
      int have = (1 << e);
      int left = (1 << e)-1;
      for (unsigned int pos = e; pos>0;) {
        int bit = 1;
        int minval = (sign ? have : -(have+left));
        int maxval = (sign ? have+left : -have);
        if (min > minval) minval = min;
        if (max < maxval) maxval = max;
        left ^= (1 << (--pos));
        if (range_test(minval,maxval)==1) return;
        int minval1 = (sign ? have+(1<<pos) : minval);
        int maxval1 = (sign ? maxval : -(have+(1<<pos)));
        int minval0 = (sign ? minval : -(have+left));
        int maxval0 = (sign ? have+left : maxval);
        if (!range_test(minval1,maxval1)) { // 1-bit is impossible
           bit = 0;
        } else if (range_test(minval0,maxval0)) { // 0-bit and 1-bit are both possible
           bit = (a >> pos) & 1;
           coder.writeMant(pos,bit);
        }
        have |= (bit << pos);
      }

    } else { // value is zero
      // only output zero bit if value could also have been nonzero
      if (range_test(min,-1) || range_test(1,max)) coder.writeZero(true);
    }
}


#endif
