OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer  -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/
CC = mpic++ -DUSE_MPI

CFLAGS = -g  $(OPTS)
CXXFLAGS = -g -x c++ --std=c++17 $(OPTS)

.PRECIOUS: %.cpt %.o

test: Makefile conjReduce
	./test.sh conjReduce

test_transformer: Makefile conjReduce.cpp ../../build/transformer
	./test_transformer.sh conjReduce.cpp

%.exe: %.tr.o Makefile inputs.o mersenne_inline.o lattice.o setup_layout_vector.o map_node_layout_trivial.o com_mpi.o memalloc.o test_gathers.o
	$(CC) -o $@ $< -lfftw3 -lm inputs.o mersenne_inline.o lattice.o setup_layout_vector.o map_node_layout_trivial.o com_mpi.o memalloc.o test_gathers.o

%.cpt: %.cpp Makefile
	$(TRANSFORMER) $<

%.cpt:../plumbing/%.cpp Makefile
	$(TRANSFORMER) $<

%.tr.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.cpt *.o *.exe
