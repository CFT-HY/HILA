
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := ../..

# Set default goal and arch
.DEFAULT_GOAL := bench

APP_HEADERS := 

APP_OPTS := -DNDIM=3

ARCH := vanilla

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

bench: build/bench ; @:

build/bench: Makefile build/bench.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/bench.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

