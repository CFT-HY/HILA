OPTS = -DNDIM=2

TRANSFORMER = ../../build/transformer
CC = clang++-8

CFLAGS = -O3 $(OPTS)
CXXFLAGS = -O3 -x c++ $(OPTS)

#add test excecutables as dependents to test and pass them to test.sh

test: Makefile conjReduce
	./test.sh conjReduce

conjReduce: Makefile conjReduce.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(CC) -o conjReduce conjReduce.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

%.cpt: %.cpp Makefile
	../../build/transformer $<

%.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf conjReduce.cpt

cleanall:
	rm -rf conjReduce.cpt conjReduce 
	rm -rf *.o
