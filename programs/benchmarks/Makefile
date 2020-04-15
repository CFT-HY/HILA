ADDITION = -DMSIZE=2
OPTS = -DNDIM=4 $(ADDITION)

#The next line work on puhti
#TRANSFORMER = ../../build/transformer -DPUHTI -DTRANSFORMER -I../../../llvm/lib/clang/8.0.1/include/ 
TRANSFORMER = ../../build/transformer  -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ 

CC = mpic++ -DUSE_MPI -mavx2 -mfma

CFLAGS = -Ofast $(OPTS)
CXXFLAGS = -Ofast -ffast-math -x c++ --std=c++17 $(OPTS)
# CXXFLAGS = -g -x c++ --std=c++17 $(OPTS)

# Don't remove cpt-files
.PRECIOUS: %.cpt

%.exe: %.tr.o Makefile inputs.o mersenne_inline.o lattice.o setup_layout_vector.o map_node_layout_trivial.o com_mpi.o timers.o memalloc.o
	$(CC) -o $@ $< inputs.o mersenne_inline.o lattice.o setup_layout_vector.o map_node_layout_trivial.o com_mpi.o timers.o memalloc.o

%.cpt: %.cpp Makefile
	$(TRANSFORMER) $(ADDITION) $<

%.tr.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@

%.tr.s : %.cpt
	$(CC) $(CXXFLAGS) $< -S -o $@

lattice.cpt : ../plumbing/lattice.cpp
	$(TRANSFORMER) $(ADDITION) $<

lattice.o : lattice.cpt
	$(CC) $(CXXFLAGS) $< -c 

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.o *.cpt *.tr.s

cleanall:
	rm -rf *.o *.cpt *.exe *.tr.s
