
# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Set default goal and arch
.DEFAULT_GOAL := hila_example

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
CUDA_VERSION=11.4
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=3

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

hila_example: build/hila_example ; @:
simple: build/simple ; @:
minmax: build/minmax; @:
# Now the linking step for each target executable
build/hila_example: Makefile build/hila_example.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/hila_example.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

# Now the linking step for each target executable
build/simple: Makefile build/simple.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/simple.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/minmax: Makefile build/minmax.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/minmax.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)






