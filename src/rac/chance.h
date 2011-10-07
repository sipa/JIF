#ifndef _RAC_CHANCE_H_
#define _RAC_CHANCE_H_ 1

#include <stdint.h>
#include <stdlib.h>

void extern build_table(uint16_t *zero_state, uint16_t *one_state, size_t size, int factor, int max_p);

class SimpleBitChanceTable {
  friend class SimpleBitChance;
private:
  int cutoff;
  uint16_t next[2][4096];
public:
  SimpleBitChanceTable(int cut) : cutoff(cut) { 
    build_table(next[0], next[1], 4096, 214748365, 4096-cut);
  }
};

extern SimpleBitChanceTable sbcTable;

class SimpleBitChance {
private:
  static const uint16_t log4k[4097];
  static const double log4k_base;

protected:
  uint16_t chance; // stored as a 12-bit number

public:
  SimpleBitChance(int chanceIn = 0x800) {
    chance = chanceIn;
  }

  uint16_t inline get() {
    return chance*16+8;
  }

  void inline put(bool bit) {
    chance = sbcTable.next[bit][chance];
  }

  void estim(bool bit, uint64_t &total) {
    total += log4k[bit ? chance : 4096-chance];
  }

  double static calcBits(uint64_t estim) {
    return log4k_base*estim;
  }
};


#endif
