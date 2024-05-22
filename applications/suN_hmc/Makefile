
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [suN_hmc]        - build suN_hmc program
#% suN_hmc make options 
#%     NCOL=<N>             - SU(N) gauge simulation program (default: 3)

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Number of colors
NCOL := 3

# Set default goal and arch
.DEFAULT_GOAL := suN_hmc

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

suN_hmc: build/suN_hmc ; @:

# Now the linking step for each target executable
build/suN_hmc: Makefile build/suN_hmc.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/suN_hmc.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 








