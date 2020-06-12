CC = mpic++ -DUSE_MPI -mavx2 -mfma
LD = mpic++

include ../../libraries/Makefile


CXXFLAGS = -Ofast -x c++ --std=c++17 $(MPI_OPTS) $(HILA_OPTS) 
#CXXFLAGS = -g -x c++ --std=c++17 $(MPI_OPTS) $(OPTS)
LDFLAGS = $(HILA_LD_FLAGS)


.PRECIOUS: build/%.cpt


staggered: Makefile build/staggered.o $(HILA_OBJECTS) staggered.h
	$(LD) $(LDFLAGS) -o build/staggered build/staggered.o $(HILA_OBJECTS)



clean:
	rm -rf build

cleanall:
	rm -rf build
