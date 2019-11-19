OPTS = -DNDIM=4

#The next line work on puhti
#TRANSFORMER = ../../build/transformer -DPUHTI -DTRANSFORMER -I../../../llvm/lib/clang/8.0.1/include/ -DNRUNS=1
TRANSFORMER = ../../build/transformer -DNRUNS=1
CC = mpic++ -DUSE_MPI -DNRUNS=1

CFLAGS = -O3 $(OPTS)
CXXFLAGS = -O3 -x c++ $(OPTS)

%.exe: %.tr.o Makefile inputs.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o  com_mpi.o
	$(CC) -o $@ $< inputs.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o

%.cpt: %.cpp Makefile
	$(TRANSFORMER) $<

%.tr.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.o *.cpt

cleanall:
	rm -rf *.o *.cpt *.exe
