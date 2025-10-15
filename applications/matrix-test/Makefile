
# Following line(s) are printed with "make help".  Use columns 8 and 30
#% matrix-test make options: (add after make [..])
#%     NCOL=<N>             - SU(N) gauge simulation program (default: 3)
#%     SUN_OVERRELAX_dFJ=1  - use deForcrand-Jahn full overrelax (default: su2 subgroups)


# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../../../HILA


# Set default goal and arch
.DEFAULT_GOAL := matrix-test

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

HILAPP_OPTS += -check-init -DEBUG

APP_OPTS += -DNDIM=3 -DBOUNDARY_LAYER_LAYOUT #-DGPU_SYNCHRONIZE_TIMERS 

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

matrix-test: build/matrix-test ; @:

# Now the linking step for each target executable
build/matrix-test: Makefile build/matrix-test.o $(HILA_OBJECTS) $(HEADERS) $(ROCM_INC_LEGACY)
	$(LD) -o $@ build/matrix-test.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) $(FFTW_LIB) -lfftw3

#build/matrix-test_bare.o: src/matrix-test_bare.h




