
TRANSFORMER = ../../build/transformer
CC = clang++-8 -g -O0 -x c++

ising: ising.cpp
	../../build/transformer -DNDIM=2 ising.cpp
	$(CC) -o ising -DNDIM=2 ising.cpt ../plumbing/lattice.cpp ../plumbing/memory.cpp ../plumbing/setup_layout_generic.cpp 