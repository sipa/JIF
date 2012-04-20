#include <vector>

#include <stdint.h>
#include "symbol.h"

// #define CONTEXT_TREE_SPLIT_THRESHOLD 5461*1
#define CONTEXT_TREE_SPLIT_THRESHOLD 5461*8*4
// k bit improvement needed before splitting

// inner nodes
class PropertyDecisionNode
{
public:
    int property;         // -1 : leaf node, childID unused
    // 0..nb_properties-1 : childID refers to left branch  (in inner_node)
    //                      childID+1 refers to right branch
    int splitval;
    int childID;
    int leafID;
    int64_t count;
    PropertyDecisionNode(int p=-1, int s=0, int c=0) : property(p), splitval(s), childID(c), leafID(0) {}
};

// leaf nodes when tree is known
template <typename BitChance> class FinalCompoundSymbolChances
{
public:
    SymbolChance<BitChance> realChances;
    int64_t count;

    void resetCounters() {
        count = 0;
    }

    FinalCompoundSymbolChances(int nBits) : realChances(nBits), count(0) {}
};

// leaf nodes during tree construction phase
template <typename BitChance> class CompoundSymbolChances : public FinalCompoundSymbolChances<BitChance>
{
public:
    std::vector<std::pair<SymbolChance<BitChance>,SymbolChance<BitChance> > > virtChances;
    uint64_t realSize;
    std::vector<uint64_t> virtSize;
    std::vector<int64_t> virtPropSum;
    signed short int best_property;

    void resetCounters() {
        FinalCompoundSymbolChances<BitChance>::resetCounters();
        best_property = -1;
        realSize = 0;
        virtPropSum.assign(virtPropSum.size(),0);
        virtSize.assign(virtSize.size(),0);
    }

    CompoundSymbolChances(int nProp, int nBits) :
        FinalCompoundSymbolChances<BitChance>(nBits),
        virtChances(nProp,std::make_pair(SymbolChance<BitChance>(nBits), SymbolChance<BitChance>(nBits))),
        realSize(0),
        virtSize(nProp),
        virtPropSum(nProp),
        best_property(-1)
    { }
};



template <typename BitChance, typename RAC> class FinalCompoundSymbolBitCoder
{
public:
    typedef typename BitChance::Table Table;

private:
    const Table &table;
    RAC &rac;
    FinalCompoundSymbolChances<BitChance> &chances;

    void inline updateChances(SymbolChanceBitType type, int i, bool bit) {
        BitChance& real = chances.realChances.bit(type,i);
        real.estim(bit, chances.realSize);
        real.put(bit, table);
    }

public:
    FinalCompoundSymbolBitCoder(const Table &tableIn, RAC &racIn, FinalCompoundSymbolChances<BitChance> &chancesIn) : table(tableIn), rac(racIn), chances(chancesIn) {}

    bool read(SymbolChanceBitType type, int i = 0) {
        BitChance& ch = chances.realChances.bit(type, i);
        bool bit = rac.read(ch.get());
        updateChances(type, i, bit);
        return bit;
    }

    void write(bool bit, SymbolChanceBitType type, int i = 0) {
        BitChance& ch = chances.realChances.bit(type, i);
        rac.write(ch.get(), bit);
        updateChances(type, i, bit);
    }
};


template <typename BitChance, typename RAC> class CompoundSymbolBitCoder
{
public:
    typedef typename BitChance::Table Table;

private:
    const Table &table;
    RAC &rac;
    CompoundSymbolChances<BitChance> &chances;
    std::vector<bool> &select;

    void inline updateChances(SymbolChanceBitType type, int i, bool bit) {
        BitChance& real = chances.realChances.bit(type,i);
        real.estim(bit, chances.realSize);
        real.put(bit, table);

        signed short int best_property = -1;
        uint64_t best_size = chances.realSize;
//    fprintf(stdout,"RealSize: %lu ||",best_size);
        for (unsigned int j=0; j<chances.virtChances.size(); j++) {
            BitChance& virt = (select)[j] ? chances.virtChances[j].first.bit(type,i)
                              : chances.virtChances[j].second.bit(type,i);
            virt.estim(bit, chances.virtSize[j]);
            virt.put(bit, table);
            if (chances.virtSize[j] < best_size) {
                best_size = chances.virtSize[j];
                best_property = j;
            }
//      fprintf(stdout,"Virt(%u)Size: %lu ||",j,chances.virtSize[j]);
        }
        chances.best_property = best_property;
//    fprintf(stdout,"\n");
    }
    BitChance inline & bestChance(SymbolChanceBitType type, int i = 0) {
        signed short int p = chances.best_property;
        return (p == -1 ? chances.realChances.bit(type,i)
                : ((select)[p] ? chances.virtChances[p].first.bit(type,i)
                   : chances.virtChances[p].second.bit(type,i) ));
    }

public:
    CompoundSymbolBitCoder(const Table &tableIn, RAC &racIn, CompoundSymbolChances<BitChance> &chancesIn, std::vector<bool> &selectIn) : table(tableIn), rac(racIn), chances(chancesIn), select(selectIn) {}

    bool read(SymbolChanceBitType type, int i = 0) {
        BitChance& ch = bestChance(type, i);
        bool bit = rac.read(ch.get());
        updateChances(type, i, bit);
//    fprintf(stderr,"bit %s%i = %s\n", SymbolChanceBitName[type], i, bit ? "true" : "false");
        return bit;
    }

    void write(bool bit, SymbolChanceBitType type, int i = 0) {
        BitChance& ch = bestChance(type, i);
        rac.write(ch.get(), bit);
        updateChances(type, i, bit);
//    fprintf(stderr,"bit %s%i = %s\n", SymbolChanceBitName[type], i, bit ? "true" : "false");
    }
};


template <typename BitChance, typename RAC> class FinalCompoundSymbolCoder
{
private:
    typedef typename FinalCompoundSymbolBitCoder<BitChance, RAC>::Table Table;
    RAC &rac;
    const Table table;

public:

    FinalCompoundSymbolCoder(RAC& racIn) : rac(racIn) {}

    int read_int(FinalCompoundSymbolChances<BitChance> &chancesIn, int min, int max) {
        FinalCompoundSymbolBitCoder<BitChance, RAC> bitCoder(table, rac, chancesIn);
        return reader(bitCoder, min, max);
    }

    void write_int(FinalCompoundSymbolChances<BitChance>& chancesIn, int min, int max, int val) {
        FinalCompoundSymbolBitCoder<BitChance, RAC> bitCoder(table, rac, chancesIn);
        writer(bitCoder, min, max, val);
    }
};

template <typename BitChance, typename RAC> class CompoundSymbolCoder
{
private:
    typedef typename CompoundSymbolBitCoder<BitChance, RAC>::Table Table;
    RAC &rac;
    const Table table;

public:

    CompoundSymbolCoder(RAC& racIn) : rac(racIn) {}

    int read_int(CompoundSymbolChances<BitChance> &chancesIn, std::vector<bool> &selectIn, int min, int max) {
        CompoundSymbolBitCoder<BitChance, RAC> bitCoder(table, rac, chancesIn, selectIn);
        return reader(bitCoder, min, max);
    }

    void write_int(CompoundSymbolChances<BitChance>& chancesIn, std::vector<bool> &selectIn, int min, int max, int val) {
        CompoundSymbolBitCoder<BitChance, RAC> bitCoder(table, rac, chancesIn, selectIn);
        writer(bitCoder, min, max, val);
    }
};



template <typename BitChance, typename RAC> class FinalPropertySymbolCoder
{
private:
    FinalCompoundSymbolCoder<BitChance, RAC> coder;
    std::vector<std::pair<int,int> > range;
    unsigned int nb_properties;
    std::vector<FinalCompoundSymbolChances<BitChance> > leaf_node;
    std::vector<PropertyDecisionNode> inner_node;

    FinalCompoundSymbolChances<BitChance> inline &find_leaf(std::vector<int> &properties) {
        int pos = 0;
        while(inner_node[pos].property != -1) {
            if (inner_node[pos].count-- > 0) break;
            if (properties[inner_node[pos].property] > inner_node[pos].splitval) {
                pos = inner_node[pos].childID;
            } else {
                pos = inner_node[pos].childID+1;
            }
        }
        FinalCompoundSymbolChances<BitChance> &result = leaf_node[inner_node[pos].leafID];

        if(inner_node[pos].count == 0) {
            int new_leaf = leaf_node.size();
            leaf_node.push_back(FinalCompoundSymbolChances<BitChance>(result));
            if (properties[inner_node[pos].property] > inner_node[pos].splitval) {
                return result;
            } else {
                return leaf_node[new_leaf];
            }
        }
        return result;
    }

public:
    FinalPropertySymbolCoder(RAC& racIn, std::vector<std::pair<int,int> > &rangeIn, int nBits, std::vector<PropertyDecisionNode> &treeIn) :
        coder(racIn),
        range(rangeIn),
        nb_properties(range.size()),
        leaf_node(1,CompoundSymbolChances<BitChance>(nb_properties,nBits)),
        inner_node(treeIn) {}

    int read_int(std::vector<int> &properties, int min, int max) {
        CompoundSymbolChances<BitChance> &chances = find_leaf(properties);
        return coder.read_int(chances, min, max);
    }

    void write_int(std::vector<int> &properties, int min, int max, int val) {
        CompoundSymbolChances<BitChance> &chances = find_leaf(properties);
        coder.write_int(chances, min, max, val);
    }

};



template <typename BitChance, typename RAC> class PropertySymbolCoder
{
public:
    typedef  std::vector<std::pair<int,int> > Ranges;
    typedef CompoundSymbolCoder<BitChance, RAC> Coder;
    typedef std::vector<PropertyDecisionNode> Tree;
private:
    RAC &rac;
    Coder coder;
    const Ranges range;
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
                pos = inner_node[pos].childID;
                //current_ranges[inner_node[pos].property].first = inner_node[pos].splitval + 1;
            } else {
                pos = inner_node[pos].childID+1;
                //current_ranges[inner_node[pos].property].second = inner_node[pos].splitval;
            }
        }
//    fprintf(stdout,"Returning leaf node %i\n", inner_node[pos].childID);
        CompoundSymbolChances<BitChance> &result = leaf_node[inner_node[pos].leafID];

        if(result.best_property != -1 && result.realSize > result.virtSize[result.best_property] + CONTEXT_TREE_SPLIT_THRESHOLD) {
            // split leaf node if some virtual context is performing (significantly) better
            int p = result.best_property;
            int new_inner = inner_node.size();
            inner_node.push_back(inner_node[pos]);
            inner_node.push_back(inner_node[pos]);
            inner_node[pos].splitval = result.virtPropSum[p]/result.count;
//            fprintf(stdout,"Splitting on property %i, splitval=%i (count=%i)\n",p,inner_node[pos].splitval, (int)result.count);
            inner_node[pos].property = p;
            inner_node[pos].count = result.count;
            int new_leaf = leaf_node.size();
            result.resetCounters();
            leaf_node.push_back(CompoundSymbolChances<BitChance>(result));
            int old_leaf = inner_node[pos].leafID;
            inner_node[pos].childID = new_inner;
            // should be OK:
            //inner_node[new_inner].childID = old_leaf;
            inner_node[new_inner].leafID = old_leaf;
            inner_node[new_inner+1].leafID = new_leaf;
            if (properties[p] > inner_node[pos].splitval) {
                return leaf_node[old_leaf];
            } else {
                return leaf_node[new_leaf];
            }
        }


        return result;
    }

    void inline set_selection_and_update_property_sums(std::vector<int> &properties, CompoundSymbolChances<BitChance> &chances) {
        chances.count++;
        for(unsigned int i=0; i<nb_properties; i++) {
            chances.virtPropSum[i] += properties[i];
//        fprintf(stdout,"Property %i: %i ||",i,properties[i]);
            selection[i] = (properties[i] > chances.virtPropSum[i]/chances.count);
        }
//    fprintf(stdout,"\n");
    }

public:
    PropertySymbolCoder(RAC& racIn, std::vector<std::pair<int,int> > &rangeIn, int nBits) :
        rac(racIn),
        coder(racIn),
        range(rangeIn),
        nb_properties(range.size()),
        leaf_node(1,CompoundSymbolChances<BitChance>(nb_properties,nBits)),
        inner_node(1,PropertyDecisionNode()),
        selection(nb_properties,false) {
    }

    int read_int(std::vector<int> &properties, int min, int max) {
        CompoundSymbolChances<BitChance> &chances = find_leaf(properties);
        set_selection_and_update_property_sums(properties,chances);
        CompoundSymbolChances<BitChance> &chances2 = find_leaf(properties);
        return coder.read_int(chances2, selection, min, max);
    }

    void write_int(std::vector<int> &properties, int min, int max, int val) {
        CompoundSymbolChances<BitChance> &chances = find_leaf(properties);
        set_selection_and_update_property_sums(properties,chances);
        CompoundSymbolChances<BitChance> &chances2 = find_leaf(properties);
        coder.write_int(chances2, selection, min, max, val);
    }

    const Tree & getTree() {
        return inner_node;
    }

};



template <typename BitChance, typename RAC> class MetaPropertySymbolCoder
{
public:
    typedef std::vector<std::pair<int,int> > Ranges;
    typedef SimpleSymbolCoder<BitChance, RAC> Coder;
    typedef std::vector<PropertyDecisionNode> Tree;
private:
    Coder coder;
    const Ranges range;
    unsigned int nb_properties;

public:
    MetaPropertySymbolCoder(RAC &racIn, const Ranges &rangesIn) :
        coder(racIn, 30),
        range(rangesIn),
        nb_properties(rangesIn.size())
    {}
    void write_subtree(int pos, Ranges &subrange, const Tree &tree) {
        PropertyDecisionNode &n = tree[pos];
        int p = n.property;
        coder.write_int(0,nb_properties,p+1);
        if (p != -1) {
            coder.write_int(0, 1<<25, n.count);
            int oldmin = subrange[p].first;
            int oldmax = subrange[p].second;
            coder.write_int(oldmin, oldmax, n.splitval);
            // > splitval
            subrange[p].first = n.splitval+1;
            write_subtree(n.childID, coder, subrange);

            // <= splitval
            subrange[p].first = oldmin;
            subrange[p].second = n.splitval;
            write_subtree(n.childID+1, coder, subrange);

            subrange[p].second = oldmax;
        }
    }
    void write_tree(const Tree &tree) {
          write_subtree(0, range, tree);
    }
    void read_subtree(int pos, Ranges &subrange, Tree &tree) {
        PropertyDecisionNode &n = tree[pos];
        int p = n.property = coder.read_int(0,nb_properties)-1;

        if (p != -1) {
            n.count = coder.read_int(0, 1<<25);
            int oldmin = subrange[p].first;
            int oldmax = subrange[p].second;
            n.splitval = coder.read_int(oldmin, oldmax);
            n.childID = tree.size();
            tree.push_back(PropertyDecisionNode());
            tree.push_back(PropertyDecisionNode());
            // > splitval
            subrange[p].first = n.splitval+1;
            read_subtree(n.childID, coder, subrange);

            // <= splitval
            subrange[p].first = oldmin;
            subrange[p].second = n.splitval;
            read_subtree(n.childID+1, coder, subrange);

            subrange[p].second = oldmax;
        }
    }
    void read_tree(Tree &tree) {
          tree.clear();
          read_subtree(0, range, tree);
    }
};