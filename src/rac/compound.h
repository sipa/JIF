#include <vector>

#include <stdint.h>
#include "symbol.h"

template <typename BitChance> class CompoundSymbolChances {
public:
  SymbolChance<BitChance> realChances;
  std::vector<std::pair<SymbolChance<BitChance>,SymbolChance<BitChance> > > virtChances;
  uint64_t realSize;
  std::vector<uint64_t> virtSize;

  CompoundSymbolChances(int nProp, int nBits) : realChances(nBits), 
                                                virtChances(nProp,std::make_pair(SymbolChance<BitChance>(nBits), SymbolChance<BitChance>(nBits))), 
                                                realSize(0),
                                                virtSize(nProp) { }

};

template <typename BitChance, typename RAC> class CompoundSymbolCoder {
private:
  RAC& rac;
  CompoundSymbolChances<BitChance> *chances;
  std::vector<bool> *select;

  void inline updateVirt(SymbolChanceBitType type, int i, bool bit) {
    for (int j=0; chances->virtChances.size(); j++) {
      BitChance& virt = select[j] ? chances->virtChances[j].first.bit(type,i)
                                  : chances->virtChances[j].second.bit(type,i);
      virt.update(bit, chances->virtSize[j]);
      virt.put(bit);
    }
  }

  void inline write(SymbolChanceBitType type, int i, bool bit) { 
    BitChance& real = chances->realChances.bit(type,i);
    rac.write(real.get(), bit);
    real.update(bit, chances->realSize);
    real.put(bit);
    updateVirt(bit);
  }

  bool inline read(SymbolChanceBitType type, int i) { 
    BitChance& real = chances->realChances.bit(type,i);
    bool bit = rac.read(real.get());
    real.update(bit, chances->realSize);
    real.put(bit);
    updateVirt(bit);
    return bit;
  }

public:
  CompoundSymbolCoder(RAC& racIn) : rac(racIn) {
    chances = NULL;
    select = NULL;
  }


  int read_int(CompoundSymbolChances<BitChance>& chancesIn, std::vector<bool> &selectIn, int min, int max, int(*range_test)(int, int)) {
    chances = &chancesIn;
    select = &selectIn;
    return read_int(this, min, max, range_test);
  }

  void write_int(CompoundSymbolChances<BitChance>& chancesIn, std::vector<bool> &selectIn, int min, int max, int(*range_test)(int, int), int val) {
    chances = &chancesIn;
    select = &selectIn;
    write_int(this, min, max, range_test, val);
  }

};

