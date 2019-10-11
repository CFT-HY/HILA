
OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer -DCUDA
CC = c++ 
#CC = nvcc -gencode arch=compute_70,code=sm_70 -DCUDA


#CUDAFLAGS = 

CFLAGS = -O3  $(OPTS) 
CXXFLAGS = -O3   $(OPTS) 
LDFLAGS = 



SUN: Makefile SUN.tr.o monte.o KennedyPendleton.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(CC) $(LDFLAGS) -o SUN SUN.tr.o monte.o KennedyPendleton.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

SUN.tr.cpp: Makefile SUN.cpp 
	$(TRANSFORMER) SUN.cpp -o=SUN.tr.cpp

SUN.tr.cu: Makefile SUN.cpp
	$(TRANSFORMER) --CUDA SUN.cpp -o=SUN.tr.cu

%.o : %.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c

%.o : %.c
	$(CC)$(CFLAGS) $< -c

clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.cpt *.o SUN

