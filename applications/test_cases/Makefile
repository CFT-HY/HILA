
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

APP_HEADERS := test.h

# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

test_CG: build/test_CG ; @:
test_conjReduce: build/test_conjReduce ; @:
test_FFT:   build/test_FFT ; @:
test_forces:   build/test_forces ; @:
bench_field:   build/bench_field ; @:
 
# Now the linking step for each target executable
build/test_CG: Makefile build/test_CG.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/test_CG.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/test_conjReduce: Makefile build/test_conjReduce.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_conjReduce.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/test_FFT: Makefile build/test_FFT.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_FFT.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/test_forces: Makefile build/test_forces.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_forces.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/bench_field: Makefile build/bench_field.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_field.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)


