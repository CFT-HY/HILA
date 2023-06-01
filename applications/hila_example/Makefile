# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Set default goal and arch
.DEFAULT_GOAL := hila_simple_example

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#Appending to hila preprocessor options
HILAPP_OPTS += -check-init

#Appending to hila applications options
APP_OPTS += -DNDIM=3

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir
hila_example: build/hila_example ; @:
hila_simple_example: build/hila_simple_example ; @:

# Now the linking step for each target executable
build/hila_example: Makefile build/hila_example.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/hila_example.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/hila_simple_example: Makefile build/hila_simple_example.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/hila_simple_example.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)