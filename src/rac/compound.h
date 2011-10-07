#include <vector>

#include <stdint.h>
#include "symbol.h"

template <typename BitChance> class CompoundSymbolChances {
public:
  SymbolChance<BitChance> realChances;
  std::vector<std::pair<SymbolChance<BitChance>,SymbolChance<BitChance> > > virtChances;
  uint64_t realSize;
  std::vector<uint64_t> virtSize;
  std::vector<uint64_t> virtPropSum;
  uint64_t count;


  CompoundSymbolChances(int nProp, int nBits) : realChances(nBits),
                                                virtChances(nProp,std::make_pair(SymbolChance<BitChance>(nBits), SymbolChance<BitChance>(nBits))),
                                                realSize(0),
                                                virtSize(nProp),
                                                virtPropSum(nProp),
                                                count(0)
                                                 { }
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
    return read_int<CompoundSymbolCoder<BitChance,RAC> >(this, min, max, range_test);
  }

  void write_int(CompoundSymbolChances<BitChance>& chancesIn, std::vector<bool> &selectIn, int min, int max, int(*range_test)(int, int), int val) {
    chances = &chancesIn;
    select = &selectIn;
    write_int<CompoundSymbolCoder<BitChance,RAC> >(this, min, max, range_test, val);
  }

};



class PropertyDecisionNode {
public:
  int property;
  int splitval;
  int left;
  int right;
  PropertyDecisionNode(int p=-1, int s=0, int l=0, int r=0) : property(p), splitval(s), left(l), right(r) {}
};

template <typename BitChance, typename RAC> class PropertySymbolCoder {
private:
  CompoundSymbolCoder<BitChance, RAC> coder;
  unsigned int nb_properties;
  std::vector<std::pair<int,int> > *range;
  std::vector<CompoundSymbolChances<BitChance> > leaf_node;
  std::vector<PropertyDecisionNode> inner_node;
  std::vector<bool> selection;

  CompoundSymbolChances<BitChance> inline &find_leaf(std::vector<int> &properties) {
    //std::vector<std:pair<int,int> > current_ranges = range;
    int pos = 0;
    while(inner_node[pos].property != -1) {
        if (properties[inner_node[pos].property] > inner_node[pos].splitval) {
                pos = inner_node[pos].left;
                //current_ranges[inner_node[pos].property].first = inner_node[pos].splitval + 1;
        } else {
                pos = inner_node[pos].right;
                //current_ranges[inner_node[pos].property].second = inner_node[pos].splitval;
        }
    }
    // todo: split leaf node if some virtual context is performing better
    return leaf_node[inner_node[pos].left];
  }

  void inline set_selection(std::vector<int> &properties, CompoundSymbolChances<BitChance> &chances) {
    chances.count++;
    for(int i=0; i<nb_properties;i++) {
        chances.virtPropSum[i] += properties[i];
        selection[i] = (properties[i] > chances.virtPropSum[i]/chances.count);
    }
  }

public:
  PropertySymbolCoder(RAC& racIn, std::vector<std::pair<int,int> > &rangeIn, int nBits) :
        coder(racIn),
        range(rangeIn),
        nb_properties(range->size()),
        leaf_node(1,CompoundSymbolChances<BitChance>(nb_properties,nBits)),
        inner_node(1),
        selection(nb_properties)
        { }

  int read_int(std::vector<int> &properties, int min, int max, int(*range_test)(int, int)) {
    CompoundSymbolChances<BitChance> chances = find_leaf(properties);
    set_selection(properties,chances);
    return coder.read_int(chances,selection, min, max, range_test);
  }

  void write_int(std::vector<int> &properties, int min, int max, int(*range_test)(int, int), int val) {
    CompoundSymbolChances<BitChance> chances = find_leaf(properties);
    set_selection(properties,chances);
    coder.write_int(chances,selection, min, max, range_test, val);
  }

};
