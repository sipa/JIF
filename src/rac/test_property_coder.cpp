#include <string>
#include <stdio.h>
#include <assert.h>

#include "symbol.h"
#include "rac.h"
#include "chance.h"
#include "compound.h"

std::string text = "The licenses for most software and other practical works are designed to take away your freedom to share and change the works.  By contrast, the GNU General Public License is intended to guarantee your freedom to share and change all versions of a program--to make sure it remains free software for all its users.  We, the Free Software Foundation, use the GNU General Public License for most of our software; it applies also to any other work released this way by its authors.  You can apply it to your programs, too.";

int main() {
  // encode
  {
  FILE *f = fopen("test_prop.dat","w");
  RacOutput40 rac(f);
  std::vector<std::pair<int,int> > ranges(2,std::make_pair(0,255));
  PropertySymbolCoder<SimpleBitChance, RacOutput40> coder(rac,ranges,8);

  int prev = 0;
  for (unsigned int i=0; i<text.size(); i++) {
    char curr = text[i];
    int min = 0 - prev;
    int max = 255 - prev;
    std::vector<int> properties(2,prev);
    int diff = curr - prev;
//    fprintf(stdout,"writing %i in %i..%i\n",diff,min,max);
    coder.write_int(properties, min, max, default_range_test, diff);
    prev = curr;
//    fprintf(stdout,"%c",(char) prev);
  }
  rac.flush();
  fclose(f);
  }
//  fprintf(stdout,"\nReading... \n");
  // decode
  {
  FILE *f = fopen("test_prop.dat","r");
  RacInput40 rac(f);
  std::vector<std::pair<int,int> > ranges(2,std::make_pair(0,255));
  PropertySymbolCoder<SimpleBitChance, RacInput40> coder(rac,ranges,8);

  int prev = 0;
  for (unsigned int i=0; i<text.size(); i++) {
    char curr = text[i];
    int min = 0 - prev;
    int max = 255 - prev;
    std::vector<int> properties(2,prev);
    int diff = coder.read_int(properties, min, max, default_range_test);
//    fprintf(stdout,"read %i in %i..%i\n",diff,min,max);
    prev = prev + diff;
    fprintf(stdout,"%c",(char) prev);
    if (curr != prev)
      printf("Fail: pos %i should be %i, but decoded %i\n",i,curr,prev);
    assert(curr == prev);
  }
  fclose(f);
  }

  return 0;
}



