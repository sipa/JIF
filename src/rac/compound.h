#include <vector>

#include <stdint.h>
#include "symbol.h"

#define CONTEXT_TREE_SPLIT_THRESHOLD 5461*4
// 4 bit improvement needed before splitting

template <typename BitChance> class CompoundSymbolChances {
public:
  SymbolChance<BitChance> realChances;
  std::vector<std::pair<SymbolChance<BitChance>,SymbolChance<BitChance> > > virtChances;
  uint64_t realSize;
  std::vector<uint64_t> virtSize;
  std::vector<int64_t> virtPropSum;
  int64_t count;
  signed short int best_property;

  void resetCounters() {
        count = 0;
        best_property = -1;
        realSize = 0;
        virtPropSum.assign(virtPropSum.size(),0);
  }

  CompoundSymbolChances(int nProp, int nBits) : realChances(nBits),
                                                virtChances(nProp,std::make_pair(SymbolChance<BitChance>(nBits), SymbolChance<BitChance>(nBits))),
                                                realSize(0),
                                                virtSize(nProp),
                                                virtPropSum(nProp),
                                                count(0),
                                                best_property(-1)
                                                 { }
};



template <typename BitChance, typename RAC> class CompoundSymbolCoder {
private:
  RAC& rac;
  CompoundSymbolChances<BitChance> *chances;
  std::vector<bool> *select;

  void inline updateChances(SymbolChanceBitType type, int i, bool bit) {
    BitChance& real = chances->realChances.bit(type,i);
    real.estim(bit, chances->realSize);
    real.put(bit);

    signed short int best_property = -1;
    uint64_t best_size = chances->realSize;
    for (unsigned int j=0; j<chances->virtChances.size(); j++) {
      BitChance& virt = (*select)[j] ? chances->virtChances[j].first.bit(type,i)
                                     : chances->virtChances[j].second.bit(type,i);
      virt.estim(bit, chances->virtSize[j]);
      virt.put(bit);
      if (chances->virtSize[j] < best_size) { best_size = chances->virtSize[j]; best_property = j; }
    }
    chances->best_property = best_property;
  }
  BitChance inline & bestChance(SymbolChanceBitType type, int i = 0) {
        signed short int p = chances->best_property;
        return (p == -1 ? chances->realChances.bit(type,i)
                        : ((*select)[p] ? chances->virtChances[p].first.bit(type,i)
                                        : chances->virtChances[p].second.bit(type,i) ));
  }

public:
  void inline write(bool bit, SymbolChanceBitType type, int i = 0) {
    BitChance& ch = bestChance(type, i);
    rac.write(ch.get(), bit);
//    fprintf(stdout,"Wrote %i with chance %i\n",bit,ch.get());
    updateChances(type, i, bit);
  }

  bool inline read(SymbolChanceBitType type, int i = 0) {
    BitChance& ch = bestChance(type, i);
    bool bit = rac.read(ch.get());
//    fprintf(stdout,"Read %i with chance %i\n",bit,ch.get());
    updateChances(type, i, bit);
    return bit;
  }

  CompoundSymbolCoder(RAC& racIn) : rac(racIn) {
    chances = NULL;
    select = NULL;
  }


  int read_int(CompoundSymbolChances<BitChance>& chancesIn, std::vector<bool> &selectIn, int min, int max, int(*range_test)(int, int)) {
    chances = &chancesIn;
    select = &selectIn;
    return reader(*this, min, max, range_test);
  }

  void write_int(CompoundSymbolChances<BitChance>& chancesIn, std::vector<bool> &selectIn, int min, int max, int(*range_test)(int, int), int val) {
    chances = &chancesIn;
    select = &selectIn;
    writer(*this, min, max, range_test, val);
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
  std::vector<std::pair<int,int> > range;
  unsigned int nb_properties;
  std::vector<CompoundSymbolChances<BitChance> > leaf_node;
  std::vector<PropertyDecisionNode> inner_node;
  std::vector<bool> selection;

  CompoundSymbolChances<BitChance> inline &find_leaf(std::vector<int> &properties) {
    //std::vector<std:pair<int,int> > current_ranges = range;
    int pos = 0;
    while(inner_node[pos].property != -1) {
//        fprintf(stderr,"Checking property %i (val=%i, splitval=%i)\n",inner_node[pos].property,properties[inner_node[pos].property],inner_node[pos].splitval);
        if (properties[inner_node[pos].property] > inner_node[pos].splitval) {
                pos = inner_node[pos].left;
                //current_ranges[inner_node[pos].property].first = inner_node[pos].splitval + 1;
        } else {
                pos = inner_node[pos].right;
                //current_ranges[inner_node[pos].property].second = inner_node[pos].splitval;
        }
    }
//    fprintf(stdout,"Returning leaf node %i\n", inner_node[pos].left);
    CompoundSymbolChances<BitChance> &result = leaf_node[inner_node[pos].left];

    if(result.best_property != -1 && result.realSize > result.virtSize[result.best_property] + CONTEXT_TREE_SPLIT_THRESHOLD) {
    // split leaf node if some virtual context is performing (significantly) better
        int p = result.best_property;
        inner_node[pos].splitval = result.virtPropSum[p]/result.count;
        fprintf(stdout,"Splitting on property %i, splitval=%i (count=%i)\n",p,inner_node[pos].splitval, result.count);
        inner_node[pos].property = p;
        int new_leaf = leaf_node.size();
        result.resetCounters();
        leaf_node.push_back(CompoundSymbolChances<BitChance>(result));
        inner_node[pos].right = new_leaf;
        if (properties[p] > inner_node[pos].splitval) {
                return leaf_node[inner_node[pos].left];
        } else {
                return leaf_node[new_leaf];
        }
    }


    return result;
  }

  void inline set_selection_and_update_property_sums(std::vector<int> &properties, CompoundSymbolChances<BitChance> &chances) {
    chances.count++;
    for(unsigned int i=0; i<nb_properties;i++) {
        chances.virtPropSum[i] += properties[i];
        selection[i] = (properties[i] > chances.virtPropSum[i]/chances.count);
    }
  }

public:
  PropertySymbolCoder(RAC& racIn, std::vector<std::pair<int,int> > &rangeIn, int nBits) :
        coder(racIn),
        range(rangeIn),
        nb_properties(range.size()),
        leaf_node(1,CompoundSymbolChances<BitChance>(nb_properties,nBits)),
        inner_node(1,PropertyDecisionNode()),
        selection(nb_properties,false)
        {
        }

  int read_int(std::vector<int> &properties, int min, int max, int(*range_test)(int, int)) {
    CompoundSymbolChances<BitChance> &chances = find_leaf(properties);
    set_selection_and_update_property_sums(properties,chances);
    CompoundSymbolChances<BitChance> &chances2 = find_leaf(properties);
    return coder.read_int(chances2,selection, min, max, range_test);
  }

  void write_int(std::vector<int> &properties, int min, int max, int(*range_test)(int, int), int val) {
    CompoundSymbolChances<BitChance> &chances = find_leaf(properties);
    set_selection_and_update_property_sums(properties,chances);
    CompoundSymbolChances<BitChance> &chances2 = find_leaf(properties);
    coder.write_int(chances2,selection, min, max, range_test, val);
  }

};
