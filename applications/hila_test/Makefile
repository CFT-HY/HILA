
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := ../..

# Set default goal and arch
.DEFAULT_GOAL := hila_test

APP_HEADERS := 

ARCH := vanilla

APP_OPTS := -DNDIM=3

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

hila_test: build/hila_test ; @:

# Now the linking step for each target executable
build/hila_test: Makefile build/hila_test.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/hila_test.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

