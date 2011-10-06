#include <string>
#include <stdio.h>
#include <assert.h>

#include "symbol.h"
#include "rac.h"
#include "chance.h"

std::string text = "The licenses for most software and other practical works are designed to take away your freedom to share and change the works.  By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users.  We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors.  You can apply it to your programs, too.";

int main() {
  {
  FILE *f = fopen("test.dat","w");
  RacOutput40 rac(f);
  SimpleSymbolCoder<SimpleBitChance, RacOutput40> coder(rac,8);

  int prev = 0;
  for (unsigned int i=0; i<text.size(); i++) {
    char curr = text[i];
    int diff = curr - prev;
    int min = 0 - prev;
    int max = 255 - prev;
    prev = curr;
//    fprintf(stdout,"writing %i in %i..%i\n",diff,min,max);
    write_int(coder, min, max, default_range_test, diff);
  }
  rac.flush();
  fclose(f);
  }

  fprintf(stdout,"\nReading... \n");

  {
  FILE *f = fopen("test.dat","r");
  RacInput40 rac(f);
  SimpleSymbolCoder<SimpleBitChance, RacInput40> coder(rac,8);

  int prev = 0;
  for (unsigned int i=0; i<text.size(); i++) {
    char curr = text[i];
    int min = 0 - prev;
    int max = 255 - prev;
    int diff = read_int(coder, min, max, default_range_test);
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



