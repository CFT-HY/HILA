
OPTS = 

TRANSFORMER = ../../build/transformer -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/
CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

TRANSFORMER_VECTOR =  ../../build/transformer  -DAVX -DUSE_MPI -target:AVX -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/
VECTOR_OPTS = -DAVX

TRANSFORMER_GPU = ../../build/transformer -DGPU -target:CUDA
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

CXXFLAGS = -Ofast -x c++ --std=c++17 $(OPTS)
LDFLAGS = -lfftw3 -lm


.PRECIOUS: build/%.cpt


OBJECTS = build/staggered.o build/inputs.o build/mersenne_inline.o build/lattice.o build/setup_layout_generic.o build/map_node_layout_trivial.o build/com_mpi.o build/memalloc.o build/test_gathers.o


staggered: Makefile $(OBJECTS) staggered.h
	$(LD) $(LDFLAGS) -o build/staggered  $(OBJECTS)


build/staggered.cpt: Makefile staggered.cpp
	mkdir -p build
	$(TRANSFORMER) staggered.cpp -o build/staggered.cpt


build/%.cpt: %.cpp Makefile
	$(TRANSFORMER) $< -o $@

build/%.cpt:../plumbing/%.cpp Makefile
	$(TRANSFORMER) $< -o $@

build/%.o : build/%.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@



clean:
	rm -rf build

cleanall:
	rm -rf build 
