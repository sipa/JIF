test_csymbol: rac.h chance.h chance.cpp symbol.h util.h util.c compound.h test_compound_symbol.cpp
	g++ -O3 -g0 -Wall util.c chance.cpp test_compound_symbol.cpp -o test_csymbol

test_csymbol.dbg: rac.h chance.h chance.cpp symbol.h util.h util.c compound.h test_compound_symbol.cpp
	g++ -O0 -ggdb3 -Wall util.c chance.cpp test_compound_symbol.cpp -o test_csymbol.dbg

test_symbol: rac.h chance.h chance.cpp symbol.h util.h util.c test_symbol.cpp
	g++ -O3 -g0 -Wall util.c chance.cpp test_symbol.cpp -o test_symbol

test_rac: rac.h chance.h chance.cpp symbol.h util.h util.c test_rac.cpp
	g++ -O3 -g0 -Wall util.c chance.cpp test_rac.cpp -o test_rac


test_symbol.dbg: rac.h chance.h chance.cpp symbol.h util.h util.c test_symbol.cpp
	g++ -O0 -ggdb3 -Wall util.c chance.cpp test_symbol.cpp -o test_symbol.dbg



