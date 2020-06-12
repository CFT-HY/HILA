OPTS = -DNDIM=4

HILAPATH = ../..
HILAINCLUDEPATH = $(HILAPATH)/libraries/

HILAPP = $(HILAPATH)/hilapp/build/hilapp -DUSE_MPI -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ -I$(HILAINCLUDEPATH)
CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

HILAPP_VECTOR =  $(HILAPATH)/hilapp/build/hilapp  -DAVX -DUSE_MPI -target:AVX -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/openmpi/include/ -I$(HILAINCLUDEPATH)
VECTOR_OPTS = -DNDIM=4

HILAPP_GPU = $(HILAPATH)/hilapp/build/hilapp -DGPU -target:CUDA -I$(HILAINCLUDEPATH)
NVCC = nvcc -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA 
NVLD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false


CXXFLAGS = -g -x c++ --std=c++17 $(OPTS) -I$(HILAINCLUDEPATH)
LDFLAGS = -lfftw3 -lm


<<<<<<< HEAD
%.cpt:plumbing/%.cpp Makefile
	$(HILAPP) $<
=======
.PRECIOUS: build/%.cpt build/%.o
>>>>>>> 344a4d92af8043ec62d7bfb79c03f59d58b4fe32

OBJECTS = build/inputs.o build/mersenne_inline.o build/lattice.o build/setup_layout_vector.o build/map_node_layout_trivial.o build/com_mpi.o build/memalloc.o build/test_gathers.o


build/%.exe: build/%.o Makefile $(OBJECTS)
	$(CC) -o $@ $< -lfftw3 -lm $(OBJECTS)

%.exe: build/%.o Makefile $(OBJECTS)
	$(CC) -o build/$@ $< -lfftw3 -lm $(OBJECTS)

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
