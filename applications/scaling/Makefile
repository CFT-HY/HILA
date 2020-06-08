OPTS = -DNDIM=3

HILAPATH = ../..
HILAINCLUDEPATH = $(HILAPATH)/libraries/

HILAPP = $(HILAPATH)/hilapp/build/hilapp -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ -I$(HILAINCLUDEPATH)
CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

TRANSFORMER_VECTOR =  $(HILAPATH)/hilapp/build/hilapp  -DAVX -DUSE_MPI -target:AVX -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ -I$(HILAINCLUDEPATH)
VECTOR_OPTS = -DAVX

TRANSFORMER_GPU = $(HILAPATH)/hilapp/build/hilapp -DGPU -target:CUDA -I$(HILAINCLUDEPATH)
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

CXXFLAGS = -Ofast -ffast-math -x c++ --std=c++17 $(OPTS)
#CXXFLAGS = -g -x c++ --std=c++17 $(OPTS) -I$(HILAINCLUDEPATH)
LDFLAGS = -lfftw3 -lm


.PRECIOUS: build/%.cpt

OBJECTS = build/inputs.o build/mersenne_inline.o build/lattice.o build/setup_layout_vector.o build/map_node_layout_trivial.o build/com_mpi.o build/memalloc.o build/test_gathers.o


onsim: Makefile build/onsim.o $(OBJECTS)
	$(LD) $(LDFLAGS) -o build/onsim build/onsim.o $(OBJECTS)

build/onsim.cpt: Makefile onsim.cpp
	mkdir -p build
	$(HILAPP) onsim.cpp -o build/onsim.cpt


build/%.cpt: %.cpp Makefile
	mkdir -p build
	$(HILAPP) $< -o $@

build/%.cpt: ../../libraries/plumbing/%.cpp Makefile
	mkdir -p build
	$(HILAPP) $< -o $@

build/%.o : build/%.cpt
	mkdir -p build
	$(CC) $(CXXFLAGS) $< -c -o $@


clean:
	rm -rf build/*.cpt build/*.o

cleanall:
	rm -rf build/*


output4: onsim
	mpirun -np 4 onsim.exe > output.dat
	sed -i '/[a-zA-Z]/d' output.dat
