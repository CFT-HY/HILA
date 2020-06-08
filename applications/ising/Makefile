
OPTS = -DNDIM=2

HILAPP = ../../build/hilapp
CC = g++
LD = g++


CXXFLAGS = -O3 -x c++ $(OPTS)



ising: Makefile ising.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o
	$(LD) -o ising ising.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o 

ising.cpt: Makefile ising.cpp
	../../build/hilapp ising.cpp

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

clean:
	rm -rf *.o ising.cpt

cleanall:
	rm -rf *.o ising.cpt ising

