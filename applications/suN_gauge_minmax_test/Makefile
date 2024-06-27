
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [minmax]        - build minmax program
#% minmax make options 
#%     NCOL=<N>              - SU(N) number of colors (default: 3)

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Number of colors
NCOL := 3

ifndef ARCH
ARCH := vanilla
endif

# Set default goal and arch
.DEFAULT_GOAL := minmax

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
#HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=4 -DNCOLOR=${NCOL} -DGFLOWACTION=${GFLOWS}


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

minmax: build/su${NCOL}_minmax_${ARCH} ; @:

# Now the linking step for each target executable
build/su${NCOL}_minmax_${ARCH}: Makefile build/minmax.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/minmax.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 








