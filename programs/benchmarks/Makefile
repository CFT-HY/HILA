OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer
CC = mpic++ -DUSE_MPI

CFLAGS = -O3 $(OPTS)
CXXFLAGS = -O3 -x c++ $(OPTS)

%.exe: %.tr.o Makefile mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o  com_mpi.o
	$(CC) -o $@ $< mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o

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
