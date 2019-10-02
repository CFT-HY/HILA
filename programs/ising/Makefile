
TRANSFORMER = ../../build/transformer
CC = clang++-8 -g -O0 -x c++

NDIM=2

ising: Makefile ising.cpt
	$(CC) -o ising -DNDIM=$(NDIM) ising.cpt ../plumbing/lattice.cpp ../plumbing/memory.cpp ../plumbing/setup_layout_generic.cpp 

ising.cpt: Makefile ising.cpp
	../../build/transformer -DNDIM=$(NDIM) ising.cpp
