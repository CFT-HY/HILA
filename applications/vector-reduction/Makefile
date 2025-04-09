
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [suN_gauge]       - build suN_hmc program

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ${HOME}/HILA

# Number of colors
NCOL := 4

# Set default goal and arch
.DEFAULT_GOAL := vector-reduction

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=3 -DNCOLOR=${NCOL}
APP_OBJECTS = build/multicanonical.o

#APP_HEADERS += twist_specific_methods.hpp checkpoint.h parameters.h
# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

vector-reduction: build/vector-reduction ; @:

# Now the linking step for each target executable
build/vector-reduction: Makefile build/vector-reduction.o $(APP_OBJECTS) $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/vector-reduction.o $(APP_OBJECTS) $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) $(FFTW_LIB)

#build/suN_gauge_bare.o: src/suN_gauge_bare.h
build/vector-reduction.o: src/vector-reduction.cpp 



