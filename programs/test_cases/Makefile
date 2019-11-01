OPTS = -DNDIM=3

TRANSFORMER = ../../build/transformer
CC = mpic++ -DUSE_MPI

CFLAGS = -g $(OPTS)
CXXFLAGS = -g -x c++ $(OPTS)

test: Makefile conjReduce
	./test.sh conjReduce

test_transformer: Makefile conjReduce.cpp ../../build/transformer
	./test_transformer.sh conjReduce.cpp

test_conjReduce: Makefile conjReduce.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 
	$(CC) -o conjReduce conjReduce.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o

%.exe: %.o Makefile mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o
	$(CC) -o $@ $< mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o
	
%.cpt: %.cpp Makefile
	../../build/transformer $<

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.o *.exe
