jif: maniac/*.h maniac/*.cpp image/*.h image/*.cpp transform/*.h transform/*.cpp jif.cpp jif.h jif_config.h
	g++ -std=gnu++11 -DNDEBUG -O3 -g0 -Wall maniac/util.c maniac/chance.cpp image/crc32k.c image/image.cpp image/image-png.cpp image/image-pnm.cpp image/color_range.cpp transform/yiq.cpp transform/zoom.cpp transform/factory.cpp jif.cpp -lpng -o jif

jif.prof: maniac/*.h maniac/*.cpp image/*.h image/*.cpp transform/*.h transform/*.cpp jif.cpp jif.h jif_config.h
	g++ -std=gnu++11 -DNDEBUG -O3 -g0 -pg -Wall maniac/util.c maniac/chance.cpp image/crc32k.c image/image.cpp image/image-png.cpp image/image-pnm.cpp image/color_range.cpp transform/yiq.cpp transform/zoom.cpp transform/factory.cpp jif.cpp -lpng -o jif.prof

jif.dbg: maniac/*.h maniac/*.cpp image/*.h image/*.cpp transform/*.h transform/*.cpp jif.cpp jif.h jif_config.h
	g++ -std=gnu++11 -O0 -ggdb3 -Wall maniac/util.c maniac/chance.cpp image/crc32k.c image/image.cpp image/image-png.cpp image/image-pnm.cpp image/color_range.cpp transform/yiq.cpp transform/zoom.cpp transform/factory.cpp jif.cpp -lpng -o jif.dbg
