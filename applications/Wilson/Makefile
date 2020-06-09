CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

include ../../libraries/Makefile_vector


CXXFLAGS = -Ofast -x c++ --std=c++17 $(MPI_OPTS) $(HILA_OPTS) 
#CXXFLAGS = -g -x c++ --std=c++17 $(MPI_OPTS) $(OPTS)
LDFLAGS = $(HILA_LD_FLAGS)


.PRECIOUS: build/%.cpt


Wilson: Makefile build/Wilson.o $(HILA_OBJECTS) Wilson.h
	$(LD) $(LDFLAGS) -o build/Wilson build/Wilson.o $(HILA_OBJECTS)



clean:
	rm -rf build

cleanall:
	rm -rf build

