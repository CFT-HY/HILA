OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer 
CC = mpic++ -DUSE_MPI

CFLAGS = -g $(OPTS)
CXXFLAGS = -g -x c++ $(OPTS)

test: Makefile conjReduce
	./test.sh conjReduce

test_transformer: Makefile conjReduce.cpp ../../build/transformer
	./test_transformer.sh conjReduce.cpp

test_conjReduce: Makefile inputs.o conjReduce.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o 
	$(CC) -o conjReduce inputs.o conjReduce.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o

%.exe: %.tr.o Makefile inputs.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o
	$(CC) -o $@ $< inputs.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o

%.cpt: %.cpp Makefile
	$(TRANSFORMER) $<

%.tr.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.cpt *.o *.exe
