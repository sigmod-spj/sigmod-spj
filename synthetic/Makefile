OBJS	= main.o MurmurHash3.o hasher.o loaddata.o
CFLAGS 	= -g -std=gnu++11

main: $(OBJS)
	g++ $(CFLAGS) $(OBJS)

main.o:	MurmurHash3.h hasher.hpp
	g++ $(CFLAGS) -c main.cpp
MurmurHash3.o: MurmurHash3.cpp
	g++ $(CFLAGS) -c MurmurHash3.cpp
hasher.o: hasher.cpp
	g++ $(CFLAGS) -c hasher.cpp
loaddata.o: loaddata.cpp
	g++ $(CFLAGS) -c loaddata.cpp
