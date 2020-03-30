
OPTS = -DNDIM=4

TRANSFORMER = ../../build/transformer 
CC = mpic++ -DUSE_MPI
LD = mpic++

TRANSFORMER_VECTOR =  ../../build/transformer -DAVX -target:AVX
VECTOR_OPTS = -DAVX

TRANSFORMER_GPU = ../../build/transformer -DGPU -target:CUDA
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

CXXFLAGS = -Ofast -ffast-math -x c++ --std=c++17 $(OPTS)
LDFLAGS = 

.PRECIOUS: %.cpt


staggered: Makefile staggered.tr.o inputs.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o
	$(LD) $(LDFLAGS) -o staggered staggered.tr.o inputs.o mersenne_inline.o lattice.o setup_layout_generic.o map_node_layout_trivial.o com_mpi.o


staggered.cpt: Makefile staggered.cpp 
	$(TRANSFORMER) staggered.cpp

%.tr.o : %.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@

%.o : %.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@
                  
%.o : ../plumbing/%.cpp
	$(CC) $(CXXFLAGS) $< -c -o $@



clean:
	rm -rf *.cpt *.o

cleanall:
	rm -rf *.cpt *.o staggered

