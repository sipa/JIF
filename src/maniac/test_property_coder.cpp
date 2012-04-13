#include <string>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <iostream>

#include "symbol.h"
#include "rac.h"
#include "chance.h"
#include "compound.h"

int main()
{
    FILE *file = fopen("test.txt", "r");
    char *text = new char[16777216];
    size_t size = fread(text, 1, 16777216, file);
    fclose(file);

    // encode
    {
        FILE *f = fopen("test_prop.dat","w");
        RacOutput40 rac(f);
        std::vector<std::pair<int,int> > ranges;
        for (int i=0; i<4; i++)
            ranges.push_back(std::make_pair(0,255));
        for (int i=0; i<3; i++)
            ranges.push_back(std::make_pair(-255,255));
        PropertySymbolCoder<SimpleBitChance, RacOutput40> coder(rac,ranges,8);

        int prev = 0;
        int prev2 = 0;
        int prev3 = 0;
        int prev4 = 0;
        for (unsigned int i=0; i<size; i++) {
            unsigned char curr = text[i];
            std::vector<int> properties(1,prev);
            properties.push_back(prev2);
            properties.push_back(prev3);
            properties.push_back(prev4);
            properties.push_back(prev - prev2);
            properties.push_back(prev2 - prev3);
            properties.push_back(prev3 - prev4);
//    fprintf(stdout,"writing %i in %i..%i\n",diff,min,max);
            coder.write_int(properties, -255, 255, default_range_test, curr - prev);
            prev4 = prev3;
            prev3 = prev2;
            prev2 = prev;
            prev = curr;
        }
        rac.flush();
        fclose(f);
    }
    printf("\n\nEncoded %li bytes\n\n", (long)size);
//  fprintf(stdout,"\nReading... \n");
    // decode
    {
        FILE *f = fopen("test_prop.dat","r");
        RacInput40 rac(f);
        std::vector<std::pair<int,int> > ranges;
        for (int i=0; i<4; i++)
            ranges.push_back(std::make_pair(0,255));
        for (int i=0; i<3; i++)
            ranges.push_back(std::make_pair(-255,255));
        PropertySymbolCoder<SimpleBitChance, RacInput40> coder(rac,ranges,8);

        int prev = 0;
        int prev2 = 0;
        int prev3 = 0;
        int prev4 = 0;
        for (unsigned int i=0; i<size; i++) {
            unsigned char curr = text[i];
            std::vector<int> properties(1,prev);
            properties.push_back(prev2);
            properties.push_back(prev3);
            properties.push_back(prev4);
            properties.push_back(prev - prev2);
            properties.push_back(prev2 - prev3);
            properties.push_back(prev3 - prev4);
            int c = coder.read_int(properties, -255, 255, default_range_test) + prev;
//    fprintf(stdout,"read %i in %i..%i\n",diff,min,max);
            prev4 = prev3;
            prev3 = prev2;
            prev2 = prev;
            prev = c;
//    fprintf(stdout,"%c",(char) c);
            if (curr != c)
                printf("Fail: pos %i should be %i, but decoded %i\n",i,curr,c);
            assert(curr == c);
        }
        fclose(f);
    }

    return 0;
}
