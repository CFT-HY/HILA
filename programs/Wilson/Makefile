
OPTS = 

TRANSFORMER = ../../build/transformer -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/
CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

TRANSFORMER_VECTOR =  ../../build/transformer  -DAVX -DUSE_MPI -target:AVX -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/
VECTOR_OPTS = -DAVX

TRANSFORMER_GPU = ../../build/transformer -DGPU -target:CUDA
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

CXXFLAGS = -Ofast -ffast-math -x c++ --std=c++17 $(OPTS)
LDFLAGS = -lfftw3 -lm


.PRECIOUS: %.cpt


Wilson: Makefile Wilson.tr.o inputs.o mersenne_inline.o lattice.o setup_layout_vector.o map_node_layout_trivial.o com_mpi.o memalloc.o test_gathers.o Wilson.h
	$(LD) $(LDFLAGS) -o Wilson Wilson.tr.o inputs.o mersenne_inline.o lattice.o setup_layout_vector.o map_node_layout_trivial.o com_mpi.o memalloc.o test_gathers.o


Wilson.cpt: Makefile Wilson.cpp 
	$(TRANSFORMER) Wilson.cpp


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
	rm -rf *.cpt *.o Wilson

