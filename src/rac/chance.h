#ifndef _RAC_CHANCE_H_
#define _RAC_CHANCE_H_ 1

#include <stdint.h>

void extern build_table(uint16_t *zero_state, uint16_t *one_state, size_t size, int factor, int max_p);

class SimpleBitChanceTable {
  friend SimpleBitChance;
private:
  int cutoff;
  uint16_t next[2][4096];
public:
  SimpleBitChanceTable(int cut) : cutoff(cut) { 
    build_table(next[0], next[1], 4096, 214748365, 4096-cut);
  }
}

extern SimpleBitChanceTable sbcTable;

class SimpleBitChance {
protected:
  uint16_t chance;

public:
  SimpleBitChance(int chanceIn = 0x8000) {
    chance = chanceIn;
  }

  int inline get() {
    return chance;
  }

  int inline put(bool bit, table_t& table) {
    chance = sbcTable.next[bit][chance];
  }
}


#endif
