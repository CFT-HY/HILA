

HILAPATH = ../..
HILAINCLUDEPATH = $(HILAPATH)/libraries/

HILAPP = $(HILAPATH)/hilapp/build/hilapp -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ -I$(HILAINCLUDEPATH)
CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

TRANSFORMER_VECTOR =  $(HILAPATH)/hilapp/build/hilapp  -DAVX -DUSE_MPI -target:AVX -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ -I$(HILAINCLUDEPATH)
VECTOR_OPTS = -DAVX

TRANSFORMER_GPU = $(HILAPATH)/hilapp/build/hilapp -DGPU -target:CUDA -I$(HILAINCLUDEPATH)
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA -I$(HILAINCLUDEPATH)
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

CXXFLAGS = -Ofast -x c++ --std=c++17 -I$(HILAINCLUDEPATH)
#CXXFLAGS = -g -x c++ --std=c++17 $(OPTS)
LDFLAGS = -lfftw3 -lm


.PRECIOUS: build/%.cpt


OBJECTS = build/Wilson.o build/inputs.o build/mersenne_inline.o build/lattice.o build/setup_layout_vector.o build/map_node_layout_trivial.o build/com_mpi.o build/memalloc.o build/test_gathers.o


Wilson: Makefile $(OBJECTS) Wilson.h
	$(LD) $(LDFLAGS) -o build/Wilson  $(OBJECTS)


build/Wilson.cpt: Makefile Wilson.cpp
	mkdir -p build
	$(HILAPP) Wilson.cpp -o build/Wilson.cpt


build/%.cpt: %.cpp Makefile
	$(HILAPP) $< -o $@

build/%.cpt: ../../libraries/plumbing/%.cpp Makefile
	$(HILAPP) $< -o $@

build/%.o : build/%.cpt
	$(CC) $(CXXFLAGS) $< -c -o $@



clean:
	rm -rf build

cleanall:
	rm -rf build

