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

static const char *SymbolChanceBitName[] = {"zero", "sign", "expo", "mant"};

static const uint16_t EXP_CHANCES[] = {3200, 2800, 2600, 2400, 2000, 1500, 800, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300};
static const uint16_t MANT_CHANCES[] = {1800, 1800, 1800, 1700, 1600, 1200, 1000, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800};
static const uint16_t ZERO_CHANCE = 1500;
static const uint16_t SIGN_CHANCE = 2048;

template <typename BitChance> class SymbolChance
{
    int bits;
    std::vector<BitChance> chances;

public:

    BitChance inline &bitZero()      {
        return chances[0];
    }
    BitChance inline &bitSign()      {
        return chances[1];
    }

    // all exp bits 0         -> int(log2(val)) == 0  [ val == 1 ]
    // exp bits up to i are 1 -> int(log2(val)) == i+1
    BitChance inline &bitExp(int i)  {
        assert(i >= 0 && i < bits-1);
        return chances[2+i];
    }
    BitChance inline &bitMant(int i) {
        assert(i >= 0 && i < bits);
        return chances[bits+1+i];
    }

    BitChance inline &bit(SymbolChanceBitType typ, int i = 0) {
        switch (typ) {
        default:
        case BIT_ZERO:
            return bitZero();
        case BIT_SIGN:
            return bitSign();
        case BIT_EXP:
            return bitExp(i);
        case BIT_MANT:
            return bitMant(i);
        }
    }


    SymbolChance(int bitsin) : bits(bitsin), chances(2*bits+1) {
        bitZero().set(ZERO_CHANCE);
        bitSign().set(SIGN_CHANCE);
        for (int i=0; i<bits-1; i++) {
            bitExp(i).set(EXP_CHANCES[i]);
        }
        for (int i=0; i<bits; i++) {
            bitMant(i).set(MANT_CHANCES[i]);
        }
    }

};

template <typename SymbolCoder> int reader(SymbolCoder& coder, int min, int max)
{
    assert(min<=max);

    if (min == max) return min; // nodig?

    if (max >= 0 && min <= 0) {
        if (coder.read(BIT_ZERO)) return 0;
    }

    bool sign = true;
    if (max > 0 && min < 0) {
        sign = coder.read(BIT_SIGN);
    } else {
        if (min < 0) sign = false;
    }
    if (sign && min <= 0) min = 1;
    if (!sign && max >= 0) max = -1;

    int amin = sign ? abs(min) : abs(max);
    int amax = sign ? abs(max) : abs(min);
    int emin = ilog2(amin), emax = ilog2(amax);
//  fprintf(stderr,"min=%i,max=%i,emin:%i,emax:%i\n",min,max,emin,emax);
    int i = emin;
    for (; i < emax; i++) {
        // if exponent >i is impossible, we are done
        if ((1 << (i+1)) > amax) break;
        if (coder.read(BIT_EXP,i)) break;
    }
    int e = i;
//  fprintf(stderr,"exp=%i\n",e);
    int have = (1 << e);
    int left = (1 << e)-1;
    for (int pos = e; pos>0;) {
        int bit = 1;
        left ^= (1 << (--pos));
        int minabs1 = have | (1<<pos);
        // int maxabs1 = have | left | (1<<pos);
        // int minabs0 = have;
        int maxabs0 = have | left;
        if (minabs1 > amax) { // 1-bit is impossible
            bit = 0;
        } else if (maxabs0 >= amin) { // 0-bit and 1-bit are both possible
            bit = coder.read(BIT_MANT,pos);
        }
        have |= (bit << pos);
    }
    int a = have;
    return (sign ? a : -a);
}

template <typename SymbolCoder> void writer(SymbolCoder& coder, int min, int max, int value)
{
    assert(min<=max);
    assert(value>=min);
    assert(value<=max);

    // avoid doing anything if the value is already known (this line is optional; behavior should be identical if commented out)
    if (min == max) return;

    if (value == 0) { // value is zero
        // only output zero bit if value could also have been nonzero
        coder.write(true, BIT_ZERO);
        return;
    }

    // only output zero bit if value could also have been zero
    if (max >= 0 && min <= 0) coder.write(false,BIT_ZERO);
    int sign = (value > 0 ? 1 : 0);
    if (max > 0 && min < 0) {
        // only output sign bit if value can be both pos and neg
        if (min < 0 && max > 0) coder.write(sign,BIT_SIGN);
    }
    if (sign && min <= 0) min = 1;
    if (!sign && max >= 0) max = -1;
    const int a = abs(value);
    const int e = ilog2(a);
    int amin = sign ? abs(min) : abs(max);
    int amax = sign ? abs(max) : abs(min);
    int emin = ilog2(amin), emax = ilog2(amax);
    // fprintf(stderr,"amin=%i,amax=%i,emin:%i,emax:%i\n",min,max,emin,emax);
    int i = emin;
    while (i < emax) {
        // if exponent >i is impossible, we are done
        if ((1 << (i+1)) > amax) break;
        // if exponent i is possible, output the exponent bit
        coder.write(i==e, BIT_EXP, i);
        if (i==e) break;
        i++;
    }
//  fprintf(stderr,"exp=%i\n",e);
    int have = (1 << e);
    int left = (1 << e)-1;
    for (int pos = e; pos>0;) {
        int bit = 1;
        left ^= (1 << (--pos));
        int minabs1 = have | (1<<pos);
        // int maxabs1 = have | left | (1<<pos);
        // int minabs0 = have;
        int maxabs0 = have | left;
        if (minabs1 > amax) { // 1-bit is impossible
            bit = 0;
        } else if (maxabs0 >= amin) { // 0-bit and 1-bit are both possible
            bit = (a >> pos) & 1;
            coder.write(bit, BIT_MANT, pos);
        }
        have |= (bit << pos);
    }
}

template <typename BitChance, typename RAC> class SimpleSymbolBitCoder
{
    typedef typename BitChance::Table Table;

private:
    const Table &table;
    SymbolChance<BitChance> &ctx;
    RAC &rac;

public:
    SimpleSymbolBitCoder(const Table &tableIn, SymbolChance<BitChance> &ctxIn, RAC &racIn) : table(tableIn), ctx(ctxIn), rac(racIn) {}

    void write(bool bit, SymbolChanceBitType typ, int i = 0) {
        BitChance& bch = ctx.bit(typ,i);
        rac.write(bch.get(), bit);
        bch.put(bit, table);
//    fprintf(stderr,"bit %s%i = %s\n", SymbolChanceBitName[typ], i, bit ? "true" : "false");
    }

    bool read(SymbolChanceBitType typ, int i = 0) {
        BitChance& bch = ctx.bit(typ,i);
        bool bit = rac.read(bch.get());
        bch.put(bit, table);
//    fprintf(stderr,"bit %s%i = %s\n", SymbolChanceBitName[typ], i, bit ? "true" : "false");
        return bit;
    }
};

template <typename BitChance, typename RAC> class SimpleSymbolCoder
{
    typedef typename BitChance::Table Table;

private:
    SymbolChance<BitChance> ctx;
    const Table table;
    RAC &rac;

public:
    SimpleSymbolCoder(RAC& racIn, int nBits) : ctx(nBits), rac(racIn) {}

    void write_int(int min, int max, int value) {
        SimpleSymbolBitCoder<BitChance, RAC> bitCoder(table, ctx, rac);
        writer(bitCoder, min, max, value);
    }

    int read_int(int min, int max) {
        SimpleSymbolBitCoder<BitChance, RAC> bitCoder(table, ctx, rac);
        return reader(bitCoder, min, max);
    }
};

#endif
