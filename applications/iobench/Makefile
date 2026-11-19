#%     make [..] NDIM=<d>    - Dimensions of system, default 4

# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := ../..

# Set default goal and arch
.DEFAULT_GOAL := iobench

APP_HEADERS :=

ARCH := vanilla

NDIM := 4

APP_OPTS := -DNDIM=$(NDIM)

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

iobench: build/iobench ; @:

build/iobench: Makefile build/iobench.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/iobench.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
