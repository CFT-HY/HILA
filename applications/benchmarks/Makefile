# Define compiler
CC = mpic++
LD = mpic++

# Define compilation flags
CXXFLAGS = -O3 -x c++ --std=c++17 
#CXXFLAGS = -g -x c++ --std=c++17 
LDFLAGS = $(HILA_LD_FLAGS)


.PRECIOUS: build/%.cpt

# Include one of the Hila library makefiles. These define the
# compilation instructions for Hila source files and define 
# $(HILA_OBJECTS), a list of objects to include in linking

# Use this one for AVX vectorization
#include ../../libraries/Makefile_vector

# Use this one for unvectorized code
include ../../libraries/Makefile

# Now the linking step for each target executable
bench_fermion: Makefile build/bench_fermion.o $(HILA_OBJECTS) bench.h
	$(LD) $(LDFLAGS) -o build/bench_fermion build/bench_fermion.o $(HILA_OBJECTS)

bench_matrix2: Makefile build/bench_matrix2.o $(HILA_OBJECTS) bench.h
	$(LD) $(LDFLAGS) -o build/bench_matrix2 build/bench_matrix2.o $(HILA_OBJECTS)

bench_field: Makefile build/bench_field.o $(HILA_OBJECTS) bench.h
	$(LD) $(LDFLAGS) -o build/bench_field build/bench_field.o $(HILA_OBJECTS)


clean:
	rm -rf build

cleanall:
	rm -rf build

