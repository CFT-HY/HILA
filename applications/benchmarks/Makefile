
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := ../..

ifndef ARCH
  $(info Defaulting to ARCH=vanilla)
  ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

APP_HEADERS := bench.h


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

bench_fermion: build/bench_fermion ; @:
bench_matrix2: build/bench_matrix2 ; @:
bench_field:   build/bench_field ; @:
bench_FFT:   build/bench_FFT ; @:

# Now the linking step for each target executable
build/bench_fermion: Makefile build/bench_fermion.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_fermion.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/bench_matrix2: Makefile build/bench_matrix2.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_matrix2.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/bench_field: Makefile build/bench_field.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_field.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/bench_FFT: Makefile build/bench_FFT.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_FFT.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)



