
OPTS = -DNDIM=2

TRANSFORMER = ../../build/transformer
CC = clang++-8

CFLAGS = -g -O0 $(OPTS)
CXXFLAGS = -g -O0 -x c++ $(OPTS)



ising: Makefile ising.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(CC) -o ising ising.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

ising.cpt: Makefile ising.cpp
	../../build/transformer ising.cpp

mersenne_inline.o: mersenne_inline.c
	gcc -c mersenne_inline.c

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf ising.cpt

cleanall:
	rm -rf ising.cpt ising

