
# Following line(s) are printed with "make help".  Use columns 8 and 30
#% suN_gauge make options:
#%     make [..] NCOL=<N>  - SU(N) gauge simulation program (default: 3)


# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Number of colors
NCOL := 3

# Set default goal and arch
.DEFAULT_GOAL := suN_gauge

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk


#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=4 -DNCOLOR=${NCOL}

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

suN_gauge: build/suN_gauge ; @:

# Now the linking step for each target executable
build/suN_gauge: Makefile build/suN_gauge.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/suN_gauge.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) $(FFTW_LIB) -lfftw3

#build/suN_gauge_bare.o: src/suN_gauge_bare.h




