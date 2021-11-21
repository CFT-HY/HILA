# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := ../..

APP_OPTS := -DNDIM=3
#-DEVEN_SITES_FIRST=0

# Set default goal and arch
.DEFAULT_GOAL := onsim

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

onsim: build/onsim ; @:

# Now the linking step for each target executable
build/onsim: Makefile build/onsim.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/onsim.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

