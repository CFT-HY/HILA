CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

OPTS = -DNDIM=3

include ../../libraries/Makefile_vector

CXXFLAGS = -O3 -x c++ --std=c++17 $(MPI_OPTS) $(HILA_OPTS) 
#CXXFLAGS = -g -x c++ --std=c++17 $(MPI_OPTS) $(OPTS)
LDFLAGS = $(HILA_LD_FLAGS)


.PRECIOUS: build/%.cpt


onsim: Makefile build/onsim.o $(HILA_OBJECTS)
	$(LD) $(LDFLAGS) -o build/onsim build/onsim.o $(HILA_OBJECTS)



clean:
	rm -rf build

cleanall:
	rm -rf build
