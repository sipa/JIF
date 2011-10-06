#include <string.h>
#include <stdio.h>

#include "symbol.h"
#include "rac.h"
#include "chance.h"

char* text = "The licenses for most software and other practical works are designed to take away your freedom to share and change the works.  By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users.  We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors.  You can apply it to your programs, too.";

int main() {
  // encode
  FILE *f = fopen("test.dat","w");
  RacOutput40 rac(f);
  SimpleSymbolCoder<SimpleBitChance, RacOutput40> coder(rac,8);

  int prev = 128; 
  for (int i=0; i<strlen(text); i++) {
    int curr = text[i];
    int diff = curr - prev;
    int min = 0 - prev;
    int max = 255 - prev;
    prev = curr;
    write_int(coder, min, max, default_range_test, diff);
  }
  rac.flush();
  fclose(f);

  // decode
  f = fopen("test.dat","r");
  RacInput40 rac(f);
  SimpleSymbolCoder<SimpleBitChance, RacInput40> decoder(rac,8);
  prev = 128;
  for (int i=0; i<strlen(text); i++) {
    int min = 0 - prev;
    int max = 255 - prev;
    int diff = read_int(coder, min, max, default_range_test);
    int curr = prev + diff;
    fprintf(stdout,"%c",curr);
    prev = curr;
  }
  fclose(f);
  return 0;
}



