jif: maniac/*.h maniac/*.cpp image/*.h image/*.cpp transform/*.h transform/*.cpp jif.cpp jif.h
	g++ -DNDEBUG -O3 -g0 -Wall maniac/util.c maniac/chance.cpp image/image.cpp image/image-pnm.cpp image/color_range.cpp transform/yiq.cpp transform/factory.cpp jif.cpp -o jif

jif.dbg: maniac/*.h maniac/*.cpp image/*.h image/*.cpp transform/*.h transform/*.cpp jif.cpp jif.h
	g++ -O0 -ggdb3 -Wall maniac/util.c maniac/chance.cpp image/image.cpp image/image-pnm.cpp image/color_range.cpp transform/yiq.cpp transform/factory.cpp jif.cpp -o jif.dbg
