
OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer -DCUDA  --CUDA
#CC = c++ 
CC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -DCUDA
LD = nvcc -arch=compute_70

#CUDAFLAGS = 

CFLAGS = -O3  $(OPTS) 
CXXFLAGS = -O3   $(OPTS) 
LDFLAGS = 



SUN: Makefile SUN.tr.o monte.o KennedyPendleton.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o
	$(LD) $(LDFLAGS) -o SUN SUN.tr.o monte.o KennedyPendleton.o mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o 

SUN.tr.cpp: Makefile SUN.cpp 
	$(TRANSFORMER) SUN.cpp -o=SUN.tr.cpp

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

