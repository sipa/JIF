#include <stdio.h>
#include <assert.h>

#include "rac.h"

int main()
{
    {
        FILE *file = fopen("bla.foo","w+");
        RacOutput40 rac(file);
        rac.write(1,4,false);
        rac.write(2,4,true);
        rac.write(3,4,false);
        for (int i=1; i<16; i++) {
            rac.write(i,16,i & 1);
            rac.write((i & 2) != 0);
        }
        rac.write(true);
        rac.flush();
        fclose(file);
    }
    {
        FILE *file = fopen("bla.foo","r");
        RacInput40 rac(file);
        assert(rac.read(1,4) == false);
        assert(rac.read(2,4) == true);
        assert(rac.read(3,4) == false);
        for (int i=1; i<16; i++) {
            assert(rac.read(i,16) == (i & 1));
            assert(rac.read() == ((i & 2) != 0));
        }
        assert(rac.read() == true);
        fclose(file);
    }
    return 0;
}
