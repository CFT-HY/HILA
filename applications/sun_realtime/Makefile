
# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Set default goal and arch
.DEFAULT_GOAL := sun_realtime

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

## Dimensions and group choice SU(N)
APP_OPTS += -DNDIM=3 -DNSU=2

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

sun_realtime: build/sun_realtime ; @:

# Now the linking step for each target executable
build/sun_realtime: Makefile build/sun_realtime.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/sun_realtime.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)







