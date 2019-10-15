
OPTS = -DNDIM=2

TRANSFORMER = ../../build/transformer
CC = g++
LD = g++


CXXFLAGS = -O3 -x c++ $(OPTS)



ising: Makefile ising.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(LD) -o ising ising.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

ising.cpt: Makefile ising.cpp
	../../build/transformer ising.cpp

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

clean:
	rm -rf ising.cpt

cleanall:
	rm -rf ising.cpt ising

