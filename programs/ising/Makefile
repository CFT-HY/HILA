
TRANSFORMER = ../../build/transformer
CC = clang++-8 -g -O0

ising: ising.cpp
	../../build/transformer -DNDIM=2 ising.cpp
	$(CC) -O -x c++ -DNDIM=2 ising.cpt ../plumbing/lattice.cpp ../plumbing/memory.cpp ../plumbing/setup_layout_generic.cpp