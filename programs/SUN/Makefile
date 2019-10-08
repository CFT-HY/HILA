
OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer --vanilla-kernel
CC = clang++-8

CFLAGS = -g -O0 $(OPTS)
CXXFLAGS = -g -O0 -x c++ $(OPTS)



SUN: Makefile SUN.o monte.o KennedyPendleton.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(CC) -o SUN SUN.o monte.o KennedyPendleton.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

SUN.cpt: Makefile SUN.cpp
	$(TRANSFORMER) SUN.cpp

monte.cpt: Makefile monte.cpp
	$(TRANSFORMER) monte.cpp


mersenne_inline.o: mersenne_inline.c
	gcc -c mersenne_inline.c

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.cpt *.o SUN

