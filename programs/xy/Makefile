
OPTS = -DNDIM=2

TRANSFORMER = ../../build/transformer
CC = g++
LD = g++


CXXFLAGS = -O3 -x c++ $(OPTS)



ising: Makefile xy.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(LD) -o xy xy.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

ising.cpt: Makefile xy.cpp
	../../build/transformer xy.cpp

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

clean:
	rm -rf *.o xy.cpt

cleanall:
	rm -rf *.o xy.cpt xy

