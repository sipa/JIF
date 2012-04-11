#include <string>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include "symbol.h"
#include "rac.h"
#include "chance.h"
#include "compound.h"

int main() {
  // encode
  {
  FILE *f = fopen("test.dat","w");
  RacOutput40 rac(f);
  CompoundSymbolChances<SimpleBitChance> chances(2,8);
  std::vector<bool> selection(2, false);
  CompoundSymbolCoder<SimpleBitChance, RacOutput40> coder(rac);

  int prev = 0;
  for (unsigned int i=0; i<size; i++) {
    char curr = text[i];
    int diff = curr - prev;
    int min = 0 - prev;
    int max = 255 - prev;
    prev = curr;
//    fprintf(stdout,"writing %i in %i..%i\n",diff,min,max);
    coder.write_int(chances, selection, min, max, default_range_test, diff);
  }
  rac.flush();
  fclose(f);
  }

//  fprintf(stdout,"\nReading... \n");
  // decode
  {
  FILE *f = fopen("test.dat","r");
  RacInput40 rac(f);
  CompoundSymbolChances<SimpleBitChance> chances(2,8);
  std::vector<bool> selection(2, false);
  CompoundSymbolCoder<SimpleBitChance, RacInput40> coder(rac);

  int prev = 0;
  for (unsigned int i=0; i<text.size(); i++) {
    char curr = text[i];
    int min = 0 - prev;
    int max = 255 - prev;
    int diff = coder.read_int(chances, selection, min, max, default_range_test);
//    fprintf(stdout,"read %i in %i..%i\n",diff,min,max);
    prev = prev + diff;
    fprintf(stdout,"%c",(char) prev);
    assert(curr == prev);
    if (curr != prev)
      printf("Fail: pos %i should be %i, but decoded %i\n",i,curr,prev);
  }
  fclose(f);
  }

  return 0;
}



